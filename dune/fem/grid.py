"""
Functions for creating python modules and C++ classes for grids.

A small example - a more complete example can be found in testgrid.py

Examples:
    >>> from __future__ import print_function

    >>> # dunefempy modules
    >>> import pydunefem
    >>> print("import this module"); import grid
    import this module...

    >>> # just get the grid (only for testing - not used)
    >>> grid1 = dune.grid.leafGrid(str("../data/unitcube-1d.dgf"), str("OneDGrid") )

    >>> # get the full grid module and then the grid (module needed for grid # functions and output object)
    >>> yaspgrid2d = dune.grid.get(str("YaspGrid"), dimgrid=2)
    >>> grid2 = yaspgrid2d.LeafGrid(str("../data/unitcube-2d.dgf"))
"""

from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import sys

from types import ModuleType

from ..generator import generator

from . import space

def interpolate(grid, func, **kwargs):
    try:
        R = func.dimRange
    except:
        R = len(func)
    spaceName = kwargs.pop('space')
    mySpace=space.create(spaceName, grid, dimrange=R, **kwargs)
    return mySpace.interpolate(func, **kwargs)

class Generator(generator.Generator):
    def modifyIncludes(self, includes):
        return includes + "#include <dune/fem/gridpart/adaptiveleafgridpart.hh>\n"
    def modifyTypeName(self, typeName):
        return "Dune::Fem::AdaptiveLeafGridPart<" + typeName + ">";
myGenerator = Generator("Grid")

def getGridType(grid, **parameters):
    """Return the grid type (using a function from database.py).
    """
    return myGenerator.getTypeName(grid, **parameters)

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
            #include <dune/fem/gridpart/adaptiveleafgridpart.hh>

            #include <dune/fempy/py/grid.hh>

            typedef Dune::Fem::AdaptiveLeafGridPart<Dune::YaspGrid< 2, Dune::EquidistantCoordinates< double, 2 > >> DuneType;

            PYBIND11_PLUGIN( grid_9f32040d49116df0211fc90ec9c908db )
            {
              pybind11::module module( "grid_9f32040d49116df0211fc90ec9c908db" );
              Dune::FemPy::registerGrid< DuneType >( module );
              return module.ptr();
            }

    This would correspond to calling yaspgrid = grid.leafGrid("../data/unitcube-2d.dgf", "YaspGrid", dimgrid=2). It also binds some functions to the created module.

    Args:
        grid (string): the identifier for the grid type to use
        parameters (kwargs): parameters used for fixing the grid type

            * dimgrid=int: dimension of the grid
            * dimworld=int: dimension of the world

    Returns:
        module: the newly created grid module

    """
    module = myGenerator.getModule(grid, **parameters)
    setattr(module.LeafGrid, "_module", module)
    setattr(module.LeafGrid, "interpolate", interpolate )
    return module

def leafGrid(dgf, grid, **parameters):
    """Get a LeafGrid

    Call get() and create a C++ grid class (see grid.hh).

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

    ret = module.LeafGrid(module.readDGF(dgf))
    return ret

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
