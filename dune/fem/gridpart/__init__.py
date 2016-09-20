from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import sys
from types import ModuleType

import dune.common as common
from ...generator.generator import SimpleGenerator
from .. import space
from .. import function as gf
from dune.fem import comm

import inspect

def interpolate(grid, func, **kwargs):
    try:
        gl = len(inspect.getargspec(func)[0])
    except:
        gl = 0
    if gl == 1:   # global function
        order = kwargs.get("order",None)
        if not order:
            order = kwargs.get("polorder",None)
            if not order:
                raise ValueError("must provide order to interpolate when used with a global function")
        return interpolate(grid, grid.globalGridFunction("gf", order, func), **kwargs)
    elif gl == 2: # local function
        order = kwargs.get("order",None)
        if not order:
            order = kwargs.get("polorder",None)
            if not order:
                raise ValueError("must provide order to interpolate when used with a local function")
        return interpolate(grid, grid.localGridFunction("gf", order, func), **kwargs)
    elif gl == 0: # already a grid function
        try:
            R = func.dimRange
        except:
            R = len(func)
        spaceName = kwargs.pop('space')
        mySpace = space.create(spaceName, grid, dimrange=R, **kwargs)
        return mySpace.interpolate(func, **kwargs)

def writeVTK(grid,  name, celldata=[], pointdata=[], cellvector=[], pointvector=[],
             number=None, subsampling=None):
    if subsampling == None:
        vtk = grid.vtkWriter()
    else:
        vtk = grid.vtkWriter(subsampling)
    for df in celldata:
        df.addToVTKWriter(vtk, common.DataType.CellData)
    for df in pointdata:
        df.addToVTKWriter(vtk, common.DataType.PointData)
    for df in cellvector:
        df.addToVTKWriter(vtk, vtk.CellVectorData)
    for df in pointvector:
        df.addToVTKWriter(vtk, vtk.PointVectorData)
    if number == None:
        vtk.write(name)
    else:
        vtk.write( name, number )
    return vtk

def levelFunction(self):
    return self.localGridFunction("level", 0, gf.Levels())
def partitionFunction(self):
    return self.localGridFunction("rank", 0, gf.Partition(comm.rank))
def globalGridFunction(grid, name, order, value):
    return grid.globalGridFunction(name, order, value)
def localGridFunction(grid, name, order, value):
    return grid.localGridFunction(name, order, value)
gridFunctions = { "globalExpr": globalGridFunction,
                  "localExpr": localGridFunction}

def function(self, name, order, *args, **kwargs):
    value = None
    key = None
    for testKey in gridFunctions:
        testValue = kwargs.pop(testKey, None)
        if testValue:
            assert not value,\
                "Only one argument allowed to define grid function"
            value = testValue
            key = testKey
    assert value,\
           "Wrong parameter used to generate grid function."+\
           "Possible parameters are:\n"+\
           ", ".join( [param for param,_ in gridFunctions.items()] )
    return gridFunctions[key](self, name, order, value, *args, **kwargs)

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "interpolate", interpolate)
    setattr(cls, "writeVTK", writeVTK)
    setattr(cls, "levelFunction", levelFunction)
    setattr(cls, "partitionFunction", partitionFunction)
    setattr(cls, "function", function)

generator = SimpleGenerator("GridPart", "Dune::FemPy")

def module(includes, typeName, constructors=None, methods=None):
    includes = includes + ["dune/fempy/py/gridpart.hh"]
    moduleName = "gridpart_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, constructors, methods)
    addAttr(module, module.GridPart)
    return module


_modules = dict()

def register(**modules):
    _modules.update(modules)


def create(gridpart, *args, **kwargs):
    gridpart = importlib.import_module(_modules[gridpart])
    return gridpart.create(*args, **kwargs)

# register our own grid parts

register(Geometry = "dune.fem.gridpart.geometry")
register(Filtered = "dune.fem.gridpart.filtered")

# enable doc test

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
