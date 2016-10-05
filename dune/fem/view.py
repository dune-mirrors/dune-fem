from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib

import dune.grid.grid_generator

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("GridView", "Dune::CorePy")

def cppBool(value):
    return "true" if value else "false"


def module(includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/corepy/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]
    moduleName = "view_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
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
        assert args.__len__() and kwargs.__len__(),\
            "too many arguments passed to adaptiveLeafGridView method"
    gridModule = importlib.import_module(type(grid).__module__)
    if isinstance(grid, getattr(gridModule, "LeafGrid")):
        grid = grid.hierarchicalGrid
        gridModule = importlib.import_module(type(grid).__module__)

    if not isinstance(grid, getattr(gridModule, "HierarchicalGrid")):
        raise ValueError('Cannot only create an adaptiveLeafGridView from a DUNE grid.')

    typeName = "Dune::Fem::AdaptiveLeafGridPart< " + grid._typeName + " >::GridViewType"
    includes = grid._includes + ["dune/fem/gridpart/adaptiveleafgridpart.hh", "dune/corepy/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]

    constructor = ["[] ( " + typeName + " &self, " + grid._typeName + " &grid ) {",
                   "    Dune::FemPy::detail::addGridModificationListener( grid );",
                   "    Dune::FemPy::constructGridPart( self, grid );",
                   "  }, pybind11::keep_alive< 1, 2 >()"]

    return module(includes, typeName, [constructor]).GridView(grid)


def filteredGridView(hostGridView, contains, useFilteredIndexSet=False):
    """create a filtered grid view

    Args:
        hostGridView:        grid view to filter
        contains:            function (Element -> bool) returning whether an element is contained in the resulting grid view
        useFilteredIndexSet: build index set containing only filtered entites? (defaults to false)

    Returns:
        GridView: the constructed grid view
    """
    includes = hostGridView._includes + ["dune/fem/gridpart/filteredgridpart.hh", "dune/fem/gridpart/filter/simple.hh"]

    hostGridViewType = hostGridView._typeName
    hostGridPartType = "Dune::FemPy::GridPart< " + hostGridViewType + " >"
    filterType = "Dune::Fem::SimpleFilter< " + hostGridPartType + " >"
    typeName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >::GridViewType"

    constructor = ["[] ( " + typeName + " &self" + ", pybind11::handle hostGridView, pybind11::function contains ) {",
                   "    auto containsCpp = [ contains ] ( const " + hostGridPartType + "::Codim< 0 >::EntityType &e ) {",
                   "        return contains( e ).template cast< bool >();",
                   "      };",
                   "    " + hostGridPartType + " &hostGridPart = Dune::FemPy::gridPart< " + hostGridViewType + " >( hostGridView );",
                   "    Dune::FemPy::constructGridPart( self, hostGridPart, " + filterType + "( hostGridPart, containsCpp ) );",
                   "  }, pybind11::keep_alive< 1, 2 >()"]

    return module(includes, typeName, [constructor]).GridView(hostGridView, contains)


def geometryGridView(coordFunction):
    """convert a coordinate function into a grid view.

    Args:
        coordFunction:  coordinate function to convert

    Returns:
        GridView: the constructed grid view
    """
    includes = coordFunction._includes + ["dune/fem/gridpart/geometrygridpart.hh"]
    typeName = "Dune::Fem::GeometryGridPart< " + coordFunction._typeName + " >::GridViewType"

    constructor = ["[] ( " + typeName + " &self" + ", " + coordFunction._typeName + " &coordFunction ) {",
                   "    Dune::FemPy::constructGridPart( self, coordFunction );",
                   "  }, pybind11::keep_alive< 1, 2 >()"]

    return module(includes, typeName, [constructor]).GridView(coordFunction)


if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
