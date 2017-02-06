from __future__ import absolute_import

import ufl
import ufl.equation

# cell
# ----

def cell(dimDomain):
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



# Space
# -----

class Space(ufl.VectorElement):
    def __init__(self, dimDomain, dimRange, field="double"):
        ufl.VectorElement.__init__(self, "Lagrange", cell(dimDomain), 1, int(dimRange))
        self._field = field
    def field(self):
        return self._field

# Coefficient
# -----------

class GridCoefficient(ufl.Coefficient):
    def __init__(self, gf):
        grid = gf.grid
        dimRange = gf.dimRange
        uflSpace = Space((grid.dimGrid, grid.dimWorld), dimRange)
        ufl.Coefficient.__init__(self, uflSpace)
        self.gf = gf

class NamedCoefficient(ufl.Coefficient):
    def __init__(self, uflSpace, name):
        ufl.Coefficient.__init__(self, uflSpace)
        self.name = name
    def str(self):
        return self.name

class NamedConstant(ufl.Coefficient):
    def __init__(self, uflSpace, name):
        ufl.Constant.__init__(self, uflSpace)
        self.name = name
    def str(self):
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
