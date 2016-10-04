from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib

import dune.grid.grid_generator

from dune.generator.generator import SimpleGenerator

generator = SimpleGenerator("GridView", "Dune::CorePy")

def adaptiveLeafGridView(grid):
    module = importlib.import_module(type(grid).__module__)
    if isinstance(grid, getattr(module, "LeafGrid")):
        grid = grid.hierarchicalGrid
        module = importlib.import_module(type(grid).__module__)

    if not isinstance(grid, getattr(module, "HierarchicalGrid")):
        raise ValueError('Cannot only create an adaptiveLeafGridView from a DUNE grid.')

    gridPartName = "Dune::Fem::AdaptiveLeafGridPart< " + grid._typeName + " >";

    typeName = gridPartName + "::GridViewType"
    includes = grid._includes + ["dune/fem/gridpart/adaptiveleafgridpart.hh", "dune/corepy/grid/gridview.hh", "dune/fempy/py/grid/gridpart.hh"]
    typeHash = "view_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()

    constructor = ["[] ( " + typeName + " &self, " + grid._typeName + " &grid ) {",
                   "    Dune::FemPy::detail::addGridModificationListener( grid );",
                   "    Dune::FemPy::constructGridPart( self, grid );",
                   "  }, pybind11::keep_alive< 1, 2 >()"]

    module = generator.load(includes, typeName, typeHash, [constructor])
    dune.grid.grid_generator.addAttr(module, module.GridView)
    return module.GridView(grid)
