from __future__ import absolute_import
from functools import wraps

import numpy
import ufl
import ufl.domain
import ufl.equation

from dune.deprecate import deprecated

# cell
# ----

##########################################
### 4d patch
##########################################
def _patchufl4d():
    from ufl.sobolevspace import H1
    from ufl.finiteelement.elementlist import ufl_elements, any_cell, register_element
    from ufl.cell import num_cell_entities, cellname2facetname, _simplex_dim2cellname, _hypercube_dim2cellname

    ## check if this has been added before
    if not 'pentatope' in ufl.cell.num_cell_entities:
        # 4d-simplex
        ufl.cell.num_cell_entities["pentatope"] = (5, 10, 10, 5, 1)
        # 4d-cube
        ufl.cell.num_cell_entities["tesseract"] = (16, 32, 24, 8, 1)

        # recompute cell name to dimension mapping
        ufl.cell.cellname2dim = dict((k, len(v) - 1) for k, v in ufl.cell.num_cell_entities.items())

        ufl.cell.cellname2facetname["pentatope"] = "tetrahedron"
        ufl.cell.cellname2facetname["tesseract"] = "hexahedron"

        ufl.cell._simplex_dim2cellname[4]   = "pentatope"
        ufl.cell._hypercube_dim2cellname[4] = "tesseract"

        # add types to element lists
        ufl.finiteelement.elementlist.simplices =\
            ufl.finiteelement.elementlist.simplices + ("pentatope",)
        ufl.finiteelement.elementlist.cubes = \
            ufl.finiteelement.elementlist.cubes + ("tesseract",)
        ufl.finiteelement.elementlist.any_cell =\
                ufl.finiteelement.elementlist.any_cell + ("pentatope", "tesseract", )

        # register Lagrange again with new element type list
        ufl_elements.pop("Lagrange")
        ufl_elements.pop("CG")
        register_element("Lagrange", "CG", 0, H1, "identity", (1, None), ufl.finiteelement.elementlist.any_cell)

def cell(dimDomainOrGrid):
    if isinstance(dimDomainOrGrid,ufl.Cell):
        return dimDomainOrGrid
    try:
        dimWorld = int(dimDomainOrGrid.dimWorld)
        dimDomain = int(dimDomainOrGrid.dimGrid)
    except:
        dimDomain = dimDomainOrGrid
        if isinstance(dimDomain, tuple):
            if len(dimDomain) != 2:
                raise Exception('dimDomain tuple must contain exactly two elements.')
            dimWorld = int(dimDomain[1])
            dimDomain = dimDomain[0]
        else:
            dimWorld = int(dimDomain)
    if dimDomain == 1:
        return ufl.Cell("interval", dimWorld)
    elif dimDomain == 2:
        return ufl.Cell("triangle", dimWorld)
    elif dimDomain == 3:
        return ufl.Cell("tetrahedron", dimWorld)
    elif dimDomain == 4:
        # add 4d cell types to ufl data structures until supported by UFL
        _patchufl4d()
        return ufl.Cell("pentatope", dimWorld)
    else:
        raise NotImplementedError('UFL cell not implemented for dimension ' + str(dimDomain) + '.')


class Space(ufl.FunctionSpace):
    def __init__(self, dimDomainOrGridOrSpace, dimRange=None, field="double", scalar=False):
        if not dimRange:
            try:
                dimRange = dimDomainOrGridOrSpace.dimRange
                dimDomainOrGridOrSpace = dimDomainOrGridOrSpace.gridView
            except AttributeError:
                dimRange = 1
                scalar = True
        self.scalar = scalar
        if scalar:
            ve = ufl.FiniteElement("Lagrange", cell(dimDomainOrGridOrSpace), 1, int(dimRange))
        else:
            ve = ufl.VectorElement("Lagrange", cell(dimDomainOrGridOrSpace), 1, int(dimRange))
        domain = ufl.domain.default_domain(ve.cell())
        ufl.FunctionSpace.__init__(self,domain, ve)
        self.dimRange = dimRange
        self.field = field
        self._cell = ve.cell()
    def cell(self):
        return self._cell
    def toVectorSpace(self):
        if not self.scalar:
            return self
        else:
            return Space(self._cell,1)

class FemSpace(Space):
    def __init__(self, space, scalar=None):
        # we shouldn't get into the situation of double wrapping
        assert not isinstance(space,FemSpace)
        if scalar is None:
            self.scalar = space.scalar
        else:
            self.scalar = scalar
        Space.__init__(self,space, scalar=self.scalar)
        self.__impl__ = space
        __module__ = space.__module__
        self.FemSpaceClass = space.__class__

    def toVectorSpace(self):
        if not self.__impl__.scalar:
            return self
        else:
            return FemSpace(self.__impl__, scalar=False)
    def as_ufl(self):
        return self

    def __getattr__(self, item):
        def tocontainer(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper
        result = getattr(self.__impl__, item)
        if not isinstance(result, FemSpace) and callable(result):
            # doc = result.func.__doc__
            result = tocontainer(result)
            # result.func.__doc__ = doc
        return result
    def __repr__(self):
        return repr(self.__impl__)
    __dict__   = property(lambda self:self.__impl__.__dict__)
    __name__   = property(lambda self:self.__impl__.__name__)
    __class__  = property(lambda self:self.__impl__.__class__)

class MixedFunctionSpace(ufl.MixedFunctionSpace):
    def __init__(self, *spaces):
        ufl.MixedFunctionSpace.__init__(self, *spaces)
        self.field = self.ufl_sub_spaces()[0].field()

    def cell(self):
        return self.ufl_element().cell()

    # def field(self):
    #     return self.ufl_sub_spaces()[0].field()


def isNumber(x):
    try:
        return 0 == x*0
    except:
        return False
# the following is an adapted version of the code in fenics
class Constant(ufl.Coefficient):
    constCount = 0
    def __init__(self, value, name=None, cell=None):
        """
        Create constant-valued function with given value.

        *Arguments*
            value
                The value may be either a single scalar value, or a
                tuple/list of values for vector-valued functions, or
                nested lists or a numpy array for tensor-valued
                functions.
            cell
                Optional argument. A :py:class:`Cell
                <ufl.Cell>` which defines the geometrical
                dimensions the Constant is defined for.
            name
                Optional argument. A str which overrules the default
                name of the Constant.

        The data type Constant represents a constant value that is
        unknown at compile-time. Its values can thus be changed
        without requiring re-generation and re-compilation of C++
        code.

        *Examples of usage*

            .. code-block:: python

                p = Constant(pi/4)              # scalar
                C = Constant((0.0, -1.0, 0.0))  # constant vector

        """

        # TODO: Either take mesh instead of cell, or drop cell and let
        # grad(c) be undefined.
        if cell is not None:
            cell = ufl.as_cell(cell)
        ufl_domain = None

        array = numpy.array(value)
        rank = len(array.shape)
        floats = list(map(float, array.flat))

        # Create UFL element and initialize constant
        if rank == 0:
            ufl_element = ufl.FiniteElement("Real", cell, 0)
        elif rank == 1:
            ufl_element = ufl.VectorElement("Real", cell, 0, dim=len(floats))
        else:
            ufl_element = ufl.TensorElement("Real", cell, 0, shape=array.shape)

        # Initialize base classes
        ufl_function_space = ufl.FunctionSpace(ufl_domain, ufl_element)
        ufl.Coefficient.__init__(self, ufl_function_space)
        if name is None:
            self.name = "c"+str(Constant.constCount)
            Constant.constCount += 1
        else:
            self.name = name
        if isNumber(value):
            self._value = float(value)
        else:
            self._value = [float(v) for v in value]
        self.models = []

    def cell(self):
        return self.ufl_element().cell()

    def values(self):
        if isNumber(self._value):
            return numpy.array([self._value])
        else:
            return self._value
    def __float__(self):
        return self._value
    @property
    def value(self):
        return self._value
    @value.setter
    def value(self,v):
        self.assign(v)
    def assign(self,v):
        if isNumber(v):
            v = float(v)
        else:
            v = [float(vv) for vv in v]
        assert type(self._value) == type(v)
        self._value = v
        for m in self.models:
            if hasattr(m,self.name):
                setattr(m,self.name,v)
            else:
                m.setConstant(self,v)
    def registerModel(self,model):
        self.models += [model]
        if hasattr(model,self.name):
            setattr(model,self.name,self._value)
        else:
            model.setConstant(self,self._value)
    # def __eq__(self, other):
    #     return self._value == other._value
    # def __ne__(self, other):
    #     return not self == other
    # def __radd__(self, other):
    #     return self.__add__(other)
    # def __add__(self, other):
    #     return self._value += other._value
    # def __sub__(self, other):
    #     return self + (-other)
    # def __rsub__(self, other):
    #     return other + (-self)
    # def __rmul__(self, scalar):
    #     return self.__mul__(scalar)
    # def __mul__(self, scalar):
    #     for i in range(len(self._value)):
    #        self._value[i] *= scalar

@deprecated("replace NamedConstant with Constant - the first argument can "\
            "now also be the initial value, e.g.,a float or list/tuple of floats")
def NamedConstant(value, name=None, dimRange=None,count=None):
    if not isinstance(value,tuple) and not isinstance(value,list) and not isNumber(value):
        try:
            domainCell = value.cell()
        except AttributeError:
            domainCell = cell(value)
        if dimRange is None:
            value = 0
        else:
            value = (0,)*dimRange
    else:
        domainCell = None
    return Constant(value, cell=domainCell, name=name)

    if dimRange is None:
        constant = ufl.Constant(domainCell, count)
        constant.values = 0
    else:
        constant = ufl.VectorConstant(domainCell, dim=dimRange, count=count)
        constant.values = (0,)*dimRange
    constant.name = name
    return constant


def Parameter(domain, parameter, dimRange=None, count=None):
    if dimRange is None:
        constant = ufl.Constant(domain, count)
    else:
        constant = ufl.VectorConstant(domain, dim=dimRange, count=count)
    constant.parameter = parameter
    return constant


from ufl.indexed import Indexed
from ufl.index_combination_utils import create_slice_indices
from ufl.core.multiindex import MultiIndex
class GridIndexed(Indexed):
    def __init__(self,gc,i):
        self.scalar = True
        component = (i,)
        shape = gc.ufl_shape
        all_indices, _, _ = create_slice_indices(component, shape, gc.ufl_free_indices)
        mi = MultiIndex(all_indices)
        Indexed.__init__(self,gc,mi)
        if gc.gf.scalar:
            assert i==0
            self.__impl__ = gc.gf
            self.gf = gc
        else:
            self.__impl__ = gc.gf[i]
            self.gf = None
    def __getattr__(self, item):
        if item == "gf": return None
        result = getattr(self.__impl__, item)
        return result
    def plot(self,*args,**kwargs):
        from dune.fem.plotting import plotPointData
        plotPointData(self.__impl__,*args,**kwargs)



class GridFunction(ufl.Coefficient):
    """ This class combines a Dune grid function and a ufl Coefficient
        class. Detailed documentation can be accessed by calling
        - help(self.GridFunction)
        - help(self.Coefficient)
    """
    def __init__(self, gf, scalar=None, count=None):
        # we shouldn't get into the situation of double wrapping
        assert not isinstance(gf,GridFunction)
        try:
            gf = gf.gf
        except:
            pass
        self.gf = gf
        self.__impl__ = gf
        __module__ = self.gf.__module__
        self.GridFunctionClass = gf.__class__
        if scalar is None and gf.dimRange == 1:
            try:
                scalar = gf.scalar
            except AttributeError:
                try:
                    scalar = gf.space.scalar
                except:
                    scalar = True

        if not hasattr(gf, "gridView"):
            gf.gridView = gf.grid # this is needed to patch dune-grid gf until the 'gridView' attribute is added there as well
        uflSpace = Space((gf.gridView.dimGrid, gf.gridView.dimWorld), gf.dimRange, scalar=scalar)

        ufl.Coefficient.__init__(self, uflSpace)
    def ufl_function_space(self):
        try:
            return self.gf.space # as_ufl()
        except TypeError or AttributeError:
            return Space(self.gf.gridView,self.gf.dimRange,scalar=False)
    def toVectorCoefficient(self):
        if not self.scalar:
            return self
        else:
            return GridFunction(self.gf,scalar=False, count=-self.count())


    def as_ufl(self):
        return self
    def copy(self,name=None):
        if name is None:
            return self.gf.copy().as_ufl()
        else:
            return self.gf.copy(name).as_ufl()

    # the following methods should be implemented here to avoid using the
    # ufl versions since they can be implemented in-place.
    def __imul__(self,a):
        self.gf.mul(a)
        return self
    def __iadd__(self,other):
        self.gf.add(other)
        return self
    def __isub__(self,other):
        self.gf.sub(other)
        return self
    def __getitem__(self,i):
        if isinstance(i,int):
            return GridIndexed(self,i)
        else:
            return ufl.Coefficient.__getitem__(self,i)
    def __getattr__(self, item):
        def tocontainer(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                ret = func(*args, **kwargs)
                if item == "localFunction":
                    setattr(ret,"scalar",self.scalar)
                return ret
            return wrapper
        result = getattr(self.__impl__, item)
        if not isinstance(result, GridFunction) and callable(result):
            # doc = result.func.__doc__
            result = tocontainer(result)
            # result.func.__doc__ = doc
        return result
    def __setattr__(self, item, v):
        super(GridFunction, self).__setattr__(item, v)
        if item == 'gf' or item == '__impl__' or item == 'GridFunctionClass':
            super(GridFunction, self).__setattr__(item, v)
        else:
            propobj = getattr(self.gf.__class__, item, None)
            if isinstance(propobj, property):
                if propobj.fset is None:
                    raise AttributeError("can't set attribute")
                propobj.fset(self, v)
            else:
                super(GridFunction, self).__setattr__(item, v)
    def __repr__(self):
        return repr(self.__impl__)
    def __str__(self):
        return self.name
    __dict__   = property(lambda self:self.gf.__dict__)
    __name__   = property(lambda self:self.gf.__name__)
    __class__  = property(lambda self:self.gf.__class__)

    def __call__(self,e,x=None):
        if x is None:
            try:
                return self.__impl__(e)
            except:
                return ufl.Coefficient.__call__(self,e)
        else:
            return self.localFunction(e)(x)

    def ufl_evaluate(self, x, component, derivatives):
        assert len(derivatives) == 0 or len(derivatives) == 1 , \
                "can only evaluate up to first order derivatives of grid functions"
        if len(derivatives) == 0:
            return self.gf.localFunction(x.entity).evaluate(x.local)[component[0]]
        else:
            return self.gf.localFunction(x.entity).jacobian(x.local)[component[0]][derivatives[0]]

class DirichletBC:
    def flatten(l):
        import collections
        for el in l:
            if isinstance(el, collections.Iterable): # and not isinstance(el, basestring):
                for sub in DirichletBC.flatten(el):
                    yield sub
            else:
                yield el

    def __init__(self, functionSpace, value, subDomain=None):
        if functionSpace.scalar:
            self.functionSpace = functionSpace.toVectorSpace()
            value = ufl.as_vector( [value] )
        else:
            self.functionSpace = functionSpace
        self.value = value
        self.subDomain = subDomain
        if type(value) is list:
            self.ufl_value = value # DirichletBC.flatten(value)
            self.ufl_value = [0 if v is None else v for v in self.ufl_value]
            self.ufl_value = ufl.as_vector(self.ufl_value)
        else:
            self.ufl_value = value
        assert self.ufl_value.ufl_shape[0] == functionSpace.dimRange
    def __str__(self):
        return str(self.value)+str(self.subDomain)
    def replace(self,dictionary):
        return DirichletBC(self.functionSpace,
                  ufl.replace(self.ufl_value,dictionary),
                  self.subDomain)
class BoxDirichletBC(DirichletBC):
    def __init__(self, functionSpace, value, xL,xR, eps=1e-10):
        cond = 1
        x = ufl.SpatialCoordinate(functionSpace)
        for l,r,c in zip(xL,xR,x):
            if l is not None:
                cond *= ufl.conditional(c>l-eps,1,0)
            if r is not None:
                cond *= ufl.conditional(c<r+eps,1,0)
        DirichletBC.__init__(self,functionSpace,value,cond)

# there is an issue here that evaluating a ufl expression can
# be very slow!
# The problem seems to be that __call__ on a ufl expression
# can only return a scalar - therefore a possible grid function
# is evaluated multiple time for each component and each
# derivative - need a better implementation?
class CoordWrapper:
    def __init__(self,e,x):
        self.entity = e
        self.local = x
        self.glb = e.geometry.toGlobal(x)
    def __getitem__(self,i): return self.glb[i]

def expression2GF(grid,expression,order,name=None):
    try:
        if expression.gf is not None:
            return expression.gf
    except:
        pass
    from dune.fem.function import localFunction, uflFunction
    return uflFunction(grid, "expr" if name is None else name, order, expression)

# register markdown formatter for integrands, forms and equations to IPython

try:
    markdown = get_ipython().display_formatter.formatters['text/markdown']
    plain = get_ipython().display_formatter.formatters['text/plain']

    from .latex import equation2latex, form2latex

    markdown.for_type(ufl.Form, lambda f: "\\begin{equation*}" + form2latex(f) + "\\end{equation*}")
    markdown.for_type(ufl.equation.Equation, lambda e: "\\begin{equation*}" + equation2latex(e) + "\\end{equation*}")

    # disable warning for calling repr to UFL forms and equations
    plain.for_type(ufl.Form, lambda f, p, c: None)
    plain.for_type(ufl.equation.Equation, lambda e, p, c: None)
except Exception as e:
    pass

import ufl.geometry
from ufl.core.ufl_type import ufl_type
from ufl.classes import all_ufl_classes, terminal_classes, ufl_classes
@ufl_type()
class BoundaryId(ufl.geometry.GeometricFacetQuantity):
    """UFL boundary id: can be used in a conditional to fix the desired boundary id."""
    __slots__ = ()
    name = "facetid"
all_ufl_classes.add(BoundaryId)
