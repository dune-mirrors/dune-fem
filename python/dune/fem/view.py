from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib

import dune.grid.grid_generator

from dune.generator import Constructor, Method

def cppBool(value):
    return "true" if value else "false"


def setup(includes, typeName, *args, ctorArgs):
    postscript="""
Dune::FemPy::registerGridView ( cls );
"""
    gv = dune.grid.grid_generator.viewModule(includes+["dune/fempy/py/gridview.hh"],
                  typeName, postscript, *args).GridView(*ctorArgs)
    gv._register()
    return gv

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

    gridPartName = "Dune::Fem::AdaptiveLeafGridPart< " + grid.cppTypeName + " >"
    typeName = gridPartName # + "::GridViewType"
    includes = grid.cppIncludes + ["dune/fem/gridpart/adaptiveleafgridpart.hh"]

    # Note: AGP are constructed from the hierarchical grid like other grid
    # views so the default ctor can be used
    '''
    constructor = Constructor([grid.cppTypeName + " &grid","int"],
                 ["std::cout << 'hallo\\n';",
                  "Dune::FemPy::detail::addGridModificationListener( grid );",
                  "return Dune::FemPy::constructGridPart<"+gridPartName+">( grid );"],
                 ["pybind11::keep_alive< 1, 2 >()"])
    '''
    return setup(includes, typeName, ctorArgs=[grid])

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
    includes = hostGridView.cppIncludes + ["dune/fem/gridpart/filteredgridpart.hh", "dune/fem/gridpart/filter/simple.hh"]

    hostGridViewType = hostGridView.cppTypeName
    hostGridPartType = "Dune::FemPy::GridPart< " + hostGridViewType + " >"
    filterType = "Dune::Fem::SimpleFilter< " + hostGridPartType + " >"
    gridPartName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >"
    typeName = "Dune::Fem::FilteredGridPart< " + hostGridPartType + ", " + filterType + ", " + cppBool(useFilteredIndexSet) + " >" # ::GridViewType"
    constructor = Constructor(["pybind11::handle hostGridView", "pybind11::function contains", "int domainId"],
                              ["auto containsCpp = [ contains ] ( const " + hostGridPartType + "::Codim< 0 >::EntityType &e ) {",
                               "    return contains( e ).template cast< int >();",
                               "  };",
                               hostGridPartType + " &hostGridPart = Dune::FemPy::gridPart< " + hostGridViewType + " >( hostGridView );",
                               "return " + gridPartName + " ( hostGridPart, " + filterType + "( hostGridPart, containsCpp, domainId ) );"],
                               # "return Dune::FemPy::constructGridPart< " + gridPartName + " >( hostGridPart, " + filterType + "( hostGridPart, containsCpp, domainId ) );"],
                              ["pybind11::keep_alive< 1, 2 >()"])
    return setup(includes, typeName, constructor, ctorArgs=[hostGridView,contains,domainId])


def geometryGridView(coordFunction):
    """convert a coordinate function into a grid view.

    Args:
        coordFunction:  coordinate function to convert

    Returns:
        GridView: the constructed grid view
    """
    assert not coordFunction.cppTypeName.startswith("Dune::Python::SimpleGridFunction"),\
"""at the moment the 'gridFunction' decorator does
not work with the 'geometryGridView'.
Interpolate into a discrete function space or use a
'uflFunction' if the function can be written as a ufl expression.
"""

    includes = coordFunction.cppIncludes + ["dune/fem/gridpart/geometrygridpart.hh"]
    gridPartName = "Dune::Fem::GeometryGridPart< " + coordFunction.cppTypeName + " >"
    typeName = gridPartName # + "::GridViewType"

    constructor = Constructor([coordFunction.cppTypeName + " &coordFunction"],
                 # ["return Dune::FemPy::constructGridPart<"+gridPartName+">( coordFunction );"],
                 ["return " + gridPartName + "( coordFunction );"],
                 ["pybind11::keep_alive< 1, 2 >()"])
    return setup(includes, typeName, constructor, ctorArgs=[coordFunction])
