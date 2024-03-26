from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import importlib

import dune.grid.grid_generator
from dune.grid import Partitions

from dune.generator import Constructor, Pickler

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

def adaptiveLeafGridView(grid, *args, partition=Partitions.all, **kwargs):
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

    partitionType = '' # default partition, All_Partition
    if partition == Partitions.interiorBorder:
        partitionType = ', Dune::InteriorBorder_Partition'
    else:
        assert partition == Partitions.all

    gridPartName = "Dune::Fem::AdaptiveLeafGridPart< " + grid.cppTypeName + partitionType + " >"
    typeName = gridPartName # + "::GridViewType"
    includes = grid.cppIncludes + ["dune/fem/gridpart/adaptiveleafgridpart.hh"]

    pickler = Pickler(getterBody=
      """
            auto& gv = self.cast<DuneType&>();
            std::ostringstream stream;
            Dune::Fem::StandardOutStream outStream(stream);
            gv.indexSet().write( outStream );
            pybind11::bytes s(stream.str());
            /* Return a tuple that fully encodes the state of the object */
            pybind11::dict d;
            if (pybind11::hasattr(self, "__dict__")) {
              d = self.attr("__dict__");
            }
            return pybind11::make_tuple(gv.grid(),s,d);
      """, setterBody=
      """
            if (t.size() != 3)
                throw std::runtime_error("Invalid state in AdaptGV with "+std::to_string(t.size())+"arguments!");
            pybind11::handle pyHg = t[0];
            auto& hg = pyHg.cast<typename DuneType::GridType&>();
            /* Create a new C++ instance */
            DuneType* gv = new DuneType(hg);
            pybind11::bytes state(t[1]);
            std::istringstream stream( state );
            Dune::Fem::StandardInStream inStream(stream);
            gv->indexSet().read( inStream );
            auto py_state = t[2].cast<pybind11::dict>();
            return std::make_pair(gv, py_state);
      """)
    return setup(includes, typeName, pickler, ctorArgs=[grid])

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
    pickler = Pickler(getterBody=
      """
            auto& gv = self.cast<DuneType&>();
            /* Return a tuple that fully encodes the state of the object */
            pybind11::dict d;
            if (pybind11::hasattr(self, "__dict__")) {
              d = self.attr("__dict__");
            }
            return pybind11::make_tuple(gv.gridFunction(),d);
      """, setterBody=
      """
            if (t.size() != 2)
                throw std::runtime_error("Invalid state in GeoGV with "+std::to_string(t.size())+"arguments!");
            pybind11::handle pyGF = t[0];
            const auto& gf = pyGF.cast<const typename DuneType::GridFunctionType&>();
            /* Create a new C++ instance */
            auto py_state = t[1].cast<pybind11::dict>();
            return std::make_pair(std::make_unique<DuneType>(gf), py_state);
      """)

    return setup(includes, typeName, constructor, pickler, ctorArgs=[coordFunction])
