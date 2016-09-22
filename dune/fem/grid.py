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
import inspect

from types import ModuleType

from dune.generator.generator import SimpleGenerator

# from .. import common
import dune.grid
import dune.grid.grid_generator
from dune.grid.core import *

from . import gridpart

class Generator(SimpleGenerator):
    def modifyIncludes(self, includes):
        return includes + ["dune/fem/gridpart/adaptiveleafgridpart.hh",
                           "dune/fempy/py/grid.hh"]
    def modifyTypeName(self, typeName):
        return "Dune::Fem::AdaptiveLeafGridPart<" + typeName + ">";a
    def load(self, includes, typeName, moduleName, constructors=None, methods=None):
        includes = self.modifyIncludes(includes)
        typeName = typeName.rsplit("::",1 )[0]
        typeName = self.modifyTypeName(typeName)
        return SimpleGenerator.load(self, includes, typeName, moduleName, constructors, methods)

def mark(self, marking):
    return self.hierarchicalGrid.mark(marking)
def adapt(self, *args):
    self.hierarchicalGrid.adapt(*args)
def globalRefine(self, *args):
    self.hierarchicalGrid.globalRefine(*args)
def loadBalance(self, *args):
    self.hierarchicalGrid.loadBalance(*args)
def triangulation(self):
    if self.dimGrid != 2 or self.dimWorld != 2:
        raise Exception("Grid must be 2-dimensional for use as matplotlib triangulation.")
    from matplotlib.tri import Triangulation
    x = self.coordinates()
    triangles = self.tesselate()
    return Triangulation(x[:,0], x[:,1], triangles)

    ret = LeafGrid(module.reader(constructor))
    return ret

def addAttr(module, cls):
    gridpart.addAttr(module,cls)
    setattr(cls, "mark", mark)
    setattr(cls, "adapt", adapt)
    setattr(cls, "globalRefine", globalRefine)
    setattr(cls, "loadBalance", loadBalance)
    setattr(cls, "triangulation", triangulation)

dune.grid.grid_generator.generator = Generator("Grid", "Dune::FemPy", "LeafGrid")
dune.grid.grid_generator.addAttr = addAttr
dune.grid.grid_generator.fileBase = "femgrid"

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
