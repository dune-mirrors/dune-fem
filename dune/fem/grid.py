"""
Functions for creating python modules and C++ classes for grids.

A small example - a more complete example can be found in testgrid.py

Examples:
    >>> from __future__ import print_function

    >>> # dunefempy modules
    >>> import pydunefem
    >>> print("import this module"); import grid
    import this module...
    >>> import gridfunction as gf

    >>> # just get the grid (only for testing - not used)
    >>> grid1 = dune.grid.leafGrid(str("../data/unitcube-1d.dgf"), str("OneDGrid") )

    >>> # get the full grid module and then the grid (module needed for grid # functions and output object)
    >>> yaspgrid2d = dune.grid.get(str("YaspGrid"), dimgrid=2)
    >>> grid2 = yaspgrid2d.LeafGrid(str("../data/unitcube-2d.dgf"))
"""

from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import sys

from functools import partial
from functools import update_wrapper
from types import ModuleType

from ..generator import generator
from . import gridfunction

class Generator(generator.Generator):
    def modifyIncludes(self, includes):
        return includes + "#include <dune/fem/gridpart/adaptiveleafgridpart.hh>\n"
    def modifyTypeName(self, typeName):
        return "Dune::Fem::AdaptiveLeafGridPart<" + typeName + ">";
myGenerator = Generator("Grid")

def getGlobalGridFunction(m,name,func):
    """NEW
       Define a global grid function.
    """
    funcR = { 1: m.GridFunctionExpression1,
              2: m.GridFunctionExpression2,
              3: m.GridFunctionExpression3,
              4: m.GridFunctionExpression4,
              5: m.GridFunctionExpression5,
              6: m.GridFunctionExpression6,
              7: m.GridFunctionExpression7,
              8: m.GridFunctionExpression8,
              9: m.GridFunctionExpression9 };
    dimR = func.dimR
    return funcR[dimR](name,func)
def getLocalGridFunction(m,name,func):
    """NEW
       Define a local grid function.
    """
    funcR = { 1: m.LocalGridFunctionExpression1,
              2: m.LocalGridFunctionExpression2,
              3: m.LocalGridFunctionExpression3,
              4: m.LocalGridFunctionExpression4,
              5: m.LocalGridFunctionExpression5,
              6: m.LocalGridFunctionExpression6,
              7: m.LocalGridFunctionExpression7,
              8: m.LocalGridFunctionExpression8,
              9: m.LocalGridFunctionExpression9 };
    dimR = func.dimR
    return funcR[dimR](name,func)
def getGlobal(self,name,func):
    return getGlobalGridFunction(self._module,name,func)
def getLocal(self,name,func):
    return getLocalGridFunction(self._module,name,func)
update_wrapper(getGlobal,getGlobalGridFunction)
update_wrapper(getLocal,getLocalGridFunction)

def vtkOutput(self):
  return self._module.VTKOutput(self)

#def interpolate(self, gridfunction, name, **parameters):
#    """NEW
#       Interpolates a grid function and returns it."""
#    from . import scheme
#
#    interScheme = scheme.scheme("InterpolationScheme", self, gridfunction, name, **parameters)
#    interScheme.solve()
#    return interScheme.solution()

def getGridType(grid, **parameters):
    """Return the grid type (using a function from database.py).
    """
    return myGenerator.getTypeName(grid, **parameters)

def addMethodsToGridModule(module):
    bind_getGlobal = partial(getGlobalGridFunction,module)
    bind_getLocal  = partial(getLocalGridFunction,module)
    setattr(module, "getGlobal", bind_getGlobal)
    setattr(module, "getLocal", bind_getLocal)


def get(grid, **parameters):
    """Create a grid module using the grid-database.

    This function creates a python module called grid_xxx.py where xxx is a
    number derived from the grid type.
    It does this by fetching the typedef and includes from the grid-database
    using the given arguments.

    Example
        A generated grid module::

            #include <dune/grid/yaspgrid.hh>
            #include <dune/grid/io/file/dgfparser/dgfyasp.hh>

            typedef Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > Grid;
            typedef Dune::Fem::AdaptiveLeafGridPart< Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > > > GridPart;

            BOOST_PYTHON_MODULE( gridc98ee9ab86a59240870d35a5c32ed4ab ) { registerDuneGrid(); }

    This would correspond to calling get("YaspGrid", dimgrid = 2). It also binds some functions to the created module.

    Args:
        grid (string): the identifier for the grid type to use
        parameters (kwargs): parameters used for fixing the grid type

            * dimgrid=int: dimension of the grid
            * dimworld=int: dimension of the world

    Returns:
        module: the newly created grid module

    """
    module = myGenerator.getModule(grid, **parameters)
    addMethodsToGridModule(module)
    return module

def addMethodsToLeafGrid(module,obj):
    #setattr(obj, "gridTypeName", "LeafGrid< " + module.gridPartTypeName + " >")
    setattr(obj, "_module", module)
    obj.getGlobal = getGlobal
    obj.getLocal = getLocal
    #obj.interpolate = interpolate
    obj.vtkOutput = vtkOutput

def leafGrid(dgf, grid, **parameters):
    """Get a LeafGrid

    Do get() and create a C++ grid class (see grid.hh).

    Notes:
        This is equivalent to::

            gridmodule = get(grid,parameters)
            grid = gridmodule.LeafGrid(dgf)

    Args:
        dgf (string): the dgf file to use for construction the grid
        grid (string): the identifier for the grid type to use
        parameters (kwargs): parameters used for fixing the grid type

            * dimgrid=int: dimension of the grid
            * dimworld=int: dimension of the world

    Returns:
        LeafGrid: the constructed grid

    """
    if isinstance(grid, str):
        module = get(grid, **parameters)
    elif isinstance(grid, ModuleType):
        module = grid
    else:
        raise TypeError("leafGrid: 'grid' must be either a string or a module")
    addMethodsToLeafGrid(module,module.LeafGrid)

    ret = module.LeafGrid(dgf)
    return ret

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
