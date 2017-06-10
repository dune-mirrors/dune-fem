from __future__ import absolute_import

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


class Space(ufl.VectorElement):
    def __init__(self, dimDomainOrGridOrSpace, dimRange=None, field="double"):
        if not dimRange:
            dimRange = dimDomainOrGridOrSpace.dimRange
            dimDomainOrGridOrSpace = dimDomainOrGridOrSpace.grid
        ufl.VectorElement.__init__(self, "Lagrange", cell(dimDomainOrGridOrSpace), 1, int(dimRange))
        self.dimRange = dimRange
        self._field = field

    def field(self):
        return self._field


def Coefficient(functionSpace, name=None, count=None):
    """UFL form argument type: Representation of a form coefficient."""
    coefficient = ufl.Coefficient(functionSpace, count=count)
    if name is not None:
        coefficient.name = name
    return coefficient


def Constant(domain, name=None, count=None):
    """UFL value: Represents a globally constant coefficient."""
    domain = ufl.domain.as_domain(domain)
    element = ufl.FiniteElement("Real", domain.ufl_cell(), 0)
    functionSpace = ufl.FunctionSpace(domain, element)
    return Coefficient(functionSpace, name=name, count=count)


def VectorConstant(domain, dimRange=None, name=None, count=None):
    """UFL value: Represents a globally constant vector-valued coefficient."""
    domain = ufl.domain.as_domain(domain)
    element = ufl.VectorElement("Real", domain.ufl_cell(), 0, dimRange)
    functionSpace = ufl.FunctionSpace(domain, element)
    return Coefficient(functionSpace, name=name, count=count)

class GridCoefficient(ufl.Coefficient):
    def __init__(self, gf):
        grid = gf.grid
        dimRange = gf.dimRange
        uflSpace = Space((grid.dimGrid, grid.dimWorld), dimRange)
        ufl.Coefficient.__init__(self, uflSpace)
        self.gf = gf
        self.__impl__ = gf
    @property
    def array(self):
        import numpy as np
        return np.array( self.gf, copy=False )
    def __getitem__(self,i):
        print("__getitem__",i)
        try:
          return GridCoefficient(self.gf[i])
        except:
          return ufl.Coefficient.__getitem__(self,i)
    def ufl_evaluate(self, x, component, derivatives):
        return self.gf.localFunction(x.entity).evaluate(x.xlocal)[component[0]]
    def __getattr__(self, item):
        result = getattr(self.__impl__, item)
        return result
    def __repr__(self):
        return repr(self.__impl__)
    def __str__(self):
        return self.name

class DirichletBC:
    def __init__(self, functionSpace, value, subDomain):
        self.functionSpace = functionSpace
        self.value = value
        self.subDomain = subDomain


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
