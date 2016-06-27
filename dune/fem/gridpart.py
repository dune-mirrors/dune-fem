from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import sys
from types import ModuleType

from ..generator import generator
from . import space
from . import function
from .. import femmpi

import inspect


def interpolate(grid, func, **kwargs):
    try:
        gl = len(inspect.getargspec(func)[0])
    except:
        gl = 0
    if gl == 1:   # global function
        return interpolate(grid, grid.globalGridFunction("gf", func), **kwargs)
    elif gl == 2: # local function
        return interpolate(grid, grid.localGridFunction("gf", func), **kwargs)
    elif gl == 0: # already a grid function
      try:
          R = func.dimRange
      except:
          R = len(func)
      spaceName = kwargs.pop('space')
      mySpace = space.create(spaceName, grid, dimrange=R, **kwargs)
      return mySpace.interpolate(func, **kwargs)

def writeVTK(grid,  name, celldata=[], pointdata=[], number=None):
    vtk = grid.vtkWriter()
    for df in celldata:
        df.addToVTKWriter(vtk, vtk.CellData)
    for df in pointdata:
        df.addToVTKWriter(vtk, vtk.PointData)
    if number == None:
        vtk.write(name)
    else:
        vtk.write( name, number )
    return vtk

def levelFunction(self):
    return self.localGridFunction("level", function.Levels())
def partitionFunction(self):
    return self.localGridFunction("rank", function.Partition(femmpi.comm.rank))

myGenerator = generator.Generator("GridPart")

def getGridPartType(gridpart, **parameters):
    """Return the gridpart type (using a function from database.py).
    """
    return myGenerator.getTypeName(gridpart, **parameters)

def addAttr(module, cls):
    setattr(cls, "_module", module)
    setattr(cls, "interpolate", interpolate )
    setattr(cls, "writeVTK", writeVTK )
    setattr(cls, "levelFunction", levelFunction )
    setattr(cls, "partitionFunction", partitionFunction )

def get(gp, **parameters):
    """Create a gridpart module using the gridpart-database.

    This function creates a python module called gridpart_xxx.py where xxx is a
    number derived from the gridpart type.
    It does this by fetching the typedef and includes from the gridpart-database
    using the given arguments.

    Returns:
        module: the newly created gridpart module

    """
    module = myGenerator.getModule(gp, **parameters)
    addAttr(module, module.GridPart)
    return module

def create(gp, gf, **parameters):
    """Get a GridPart.

    Call get() and create a C++ gridpart class

    Args:
        gridpart (string): the identifier for the scheme type to use
        gf; main argument, will in some form conatin the host grid part
    Returns:
        GridPart: the constructed GridPart
    """
    if gp=="Geometry":
        # try:
            module = get(gp, discfunc=gf._module._typeName, extra_includes=gf._module._includes, **parameters)
        # except:
        #     gp  = "Hallo"
        #     val = "Dune::FieldVector<double,"+str(gf.dimRange)+">"
        #     module = get("GeometryVirtual", gridpart=gp, value=val, **parameters)
    return module.GridPart(gf)

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
