from __future__ import absolute_import
from functools import wraps

import ufl
import ufl.domain
import ufl.equation

# cell
# ----

def cell(dimDomainOrGrid):
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
    else:
        raise NotImplementedError('UFL cell not implemented for dimension' + str(dimDomain) + '.')


class Space(ufl.FunctionSpace):
    def __init__(self, dimDomainOrGridOrSpace, dimRange=None, field="double"):
        if not dimRange:
            self.duneSpace = dimDomainOrGridOrSpace
            dimRange = dimDomainOrGridOrSpace.dimRange
            dimDomainOrGridOrSpace = dimDomainOrGridOrSpace.grid
        ve = ufl.VectorElement("Lagrange", cell(dimDomainOrGridOrSpace), 1, int(dimRange))
        domain = ufl.domain.default_domain(ve.cell())
        ufl.FunctionSpace.__init__(self,domain, ve)
        self.dimRange = dimRange
        self._field = field
        self._cell = ve.cell()
    def cell(self):
        return self._cell
    def field(self):
        return self._field


class MixedFunctionSpace(ufl.MixedFunctionSpace):
    def __init__(self, *spaces):
        ufl.MixedFunctionSpace.__init__(self, *spaces)

    def cell(self):
        return self.ufl_element().cell()

    def field(self):
        return self.ufl_sub_spaces()[0].field()


def NamedCoefficient(functionSpace, name=None, count=None):
    coefficient = ufl.Coefficient(functionSpace, count=count)
    if name is not None:
        coefficient.name = name
    return coefficient
def NamedConstant(domain, dimRange=None, name=None, count=None):
    if dimRange == 0:
        constant = ufl.Constant(domain, count)
    else:
        constant = ufl.VectorConstant(domain, dim=dimRange, count=count)
    if name is not None:
        constant.name = name
    return constant


from ufl.indexed import Indexed
from ufl.index_combination_utils import create_slice_indices
from ufl.core.multiindex import MultiIndex
class GridIndexed(Indexed):
    def __init__(self,gc,i):
        component = (i,)
        shape = gc.ufl_shape
        all_indices, _, _ = create_slice_indices(component, shape, gc.ufl_free_indices)
        mi = MultiIndex(all_indices)
        Indexed.__init__(self,gc,mi)
        self.__impl__ = gc.gf[i]
    def __getattr__(self, item):
        result = getattr(self.__impl__, item)
        return result



class GridFunction(ufl.Coefficient):
    """ This class combines a Dune grid function and a ufl Coefficient
        class. Detailed documentation can be accessed by calling
        - help(self.GridFunction)
        - help(self.Coefficient)
    """
    def __init__(self, gf):
        try:
            gf = gf.gf
        except:
            pass
        self.gf = gf
        self.__impl__ = gf
        __module__ = self.gf.__module__
        self.GridFunctionClass = gf.__class__
        grid = gf.grid

        dimRange = gf.dimRange
        uflSpace = Space((grid.dimGrid, grid.dimWorld), dimRange)
        ufl.Coefficient.__init__(self, uflSpace)
    def copy(self):
        return self.gf.copy().as_ufl(); # GridFunction(self.gf.copy())
    @property
    def as_numpy(self):
        import numpy as np
        return np.array( self.dofVector, copy=False )
    @property
    def array(self):
        import numpy as np
        return np.array( self.gf, copy=False )
    def __getitem__(self,i):
        if isinstance(i,int):
            return GridIndexed(self,i)
        else:
            return ufl.Coefficient.__getitem__(self,i)
    def __getattr__(self, item):
        def tocontainer(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                return func(*args, **kwargs)
            return wrapper
        result = getattr(self.__impl__, item)
        if callable(result):
            # doc = result.func.__doc__
            result = tocontainer(result)
            # result.func.__doc__ = doc
        return result
    def __repr__(self):
        return repr(self.__impl__)
    def __str__(self):
        return self.name
    __dict__   = property(lambda self:self.gf.__dict__)
    __name__   = property(lambda self:self.gf.__name__)
    __class__  = property(lambda self:self.gf.__class__)

    def ufl_evaluate(self, x, component, derivatives):
        assert len(derivatives) == 0 or len(derivatives) == 1 , \
                "can only evaluate up to first order derivatives of grid functions"
        if len(derivatives) == 0:
            return self.gf.localFunction(x.entity).evaluate(x.local)[component[0]]
        else:
            return self.gf.localFunction(x.entity).jacobian(x.local)[component[0]][derivatives[0]]
    def __getitem__(self,i):
        if isinstance(i,int):
            return GridIndexed(self,i)
        else:
            return ufl.Coefficient.__getitem__(self,i)

class DirichletBC:
    def __init__(self, functionSpace, value, subDomain):
        try:
            self.functionSpace = functionSpace.uflSpace
        except AttributeError:
            self.functionSpace = functionSpace
        self.value = value
        self.subDomain = subDomain

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
        self.glb = e.geometry.position(x)
    def __getitem__(self,i): return self.glb[i]
def expression2GF(grid,expression,order):
    from dune.fem.function._functions import localFunction
    shape = expression.ufl_shape
    assert len(shape) == 0 or len(shape) == 1,\
            "can only generate grid function from scalar or vector valued expression, got %s" %str(shape)
    if len(shape) == 0:
        return localFunction(grid, "tmp", order, lambda e,x: [expression(CoordWrapper(e,x))] )
    if len(shape) == 1:
        return localFunction(grid, "tmp", order, lambda e,x: [expression[i](CoordWrapper(e,x)) for i in range(shape[0]) ] )

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
