from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib

import dune.grid.grid_generator

from dune.generator import Constructor, Method
from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("GridView", "Dune::FemPy")

def cppBool(value):
    return "true" if value else "false"


def load(includes, typeName, *args):
    # includes = includes + ["dune/fempy/py/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]
    pyIncludes = ["dune/fempy/py/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]
    moduleName = "view_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    holder = "Dune::FemPy::GridPartPtr< " + typeName + " >"
    # module = generator.load(includes, typeName, moduleName, *args, options=[holder])
    module = generator.load([includes,pyIncludes], typeName, moduleName, *args)
    dune.grid.grid_generator.addAttr(module, module.GridView)
    return module


def adaptiveLeafGridView(grid, *args, **kwargs):
    """create an adaptive view of the leaf grid

    Args:
        grid:  grid to create the adaptive view for.
               The grid must either be a hierarchical grid or a leaf view of one.

    Returns:
        GridView: the constructed grid view
    """
    if isinstance(grid, str):
        import dune.create as create
        grid = create.grid(grid,*args,**kwargs)
    else:
        assert args.__len__()==0 and kwargs.__len__()==0,\
            "too many arguments passed to adaptiveLeafGridView method"

    try:
        grid = grid.hierarchicalGrid
    except:
        pass
    gridModule = importlib.import_module(type(grid).__module__)

    if not isinstance(grid, getattr(gridModule, "HierarchicalGrid")):
        raise ValueError('Cannot only create an adaptiveLeafGridView from a DUNE grid.')

    gridPartName = "Dune::Fem::AdaptiveLeafGridPart< " + grid._typeName + " >"
    typeName = gridPartName + "::GridViewType"
    includes = grid._includes + ["dune/fem/gridpart/adaptiveleafgridpart.hh", "dune/python/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]

    # constructor = Constructor([grid._typeName + " &grid"],
    #                           ["Dune::FemPy::detail::addGridModificationListener( grid );",
    #                            "return Dune::FemPy::makeGridPart< " + typeName + " >( grid );"],
    #                           ["pybind11::keep_alive< 1, 2 >()"])
    constructor = Constructor([grid._typeName + " &grid"],
                 ["Dune::FemPy::detail::addGridModificationListener( grid );",
                  "return Dune::FemPy::constructGridPart<"+gridPartName+">( grid );"],
                 ["pybind11::keep_alive< 1, 2 >()"])

    GridView = load(includes, typeName, constructor).GridView
    # setattr(GridView,"canAdapt",True)
    return GridView(grid)


def filteredGridView(hostGridView, contains, domainId, useFilteredIndexSet=False):
    """create a filtered grid view

    Args:
        hostGridView:        grid view to filter
        contains:            function (Element -> int) returns a domain id for each element is contained in the resulting grid view
        domainId:            contains==domainId used to define entities inside the filtered gv
        useFilteredIndexSet: build index set containing only filtered entites? (defaults to false)

    Returns:
        GridView: the constructed grid view
    """
    includes = hostGridView._includes + ["dune/fem/gridpart/filteredgridpart.hh", "dune/fem/gridpart/filter/simple.hh", "dune/python/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]

    hostGridViewType = hostGridView._typeName
    hostGridPartType = "Dune::FemPy::GridPart< " + hostGridViewType + " >"
    filterType = "Dune::Fem::SimpleFilter< " + hostGridPartType + " >"
    gridPartName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >"
    typeName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >::GridViewType"

    constructor = Constructor(["pybind11::handle hostGridView", "pybind11::function contains", "int domainId"],
                              ["auto containsCpp = [ contains ] ( const " + hostGridPartType + "::Codim< 0 >::EntityType &e ) {",
                               "    return contains( e ).template cast< int >();",
                               "  };",
                               hostGridPartType + " &hostGridPart = Dune::FemPy::gridPart< " + hostGridViewType + " >( hostGridView );",
                               "return Dune::FemPy::constructGridPart< " + gridPartName + " >( hostGridPart, " + filterType + "( hostGridPart, containsCpp, domainId ) );"],
                              ["pybind11::keep_alive< 1, 2 >()"])

    return load(includes, typeName, constructor).GridView(hostGridView, contains, domainId)


def geometryGridView(coordFunction):
    """convert a coordinate function into a grid view.

    Args:
        coordFunction:  coordinate function to convert

    Returns:
        GridView: the constructed grid view
    """
    includes = coordFunction._includes + ["dune/fem/gridpart/geometrygridpart.hh", "dune/python/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]
    gridPartName = "Dune::Fem::GeometryGridPart< " + coordFunction._typeName + " >"
    typeName = gridPartName + "::GridViewType"

    # constructor = Constructor([coordFunction._typeName + " &coordFunction"],
    #                           ["return Dune::FemPy::makeGridPart< " + typeName + " >( coordFunction );"],
    #                           ["pybind11::keep_alive< 1, 2 >()"])
    constructor = Constructor([coordFunction._typeName + " &coordFunction"],
                 ["return Dune::FemPy::constructGridPart<"+gridPartName+">( coordFunction );"],
                 ["pybind11::keep_alive< 1, 2 >()"])

    return load(includes, typeName, constructor).GridView(coordFunction)


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
