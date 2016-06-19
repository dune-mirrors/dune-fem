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

myGenerator = generator.Generator("GridPart")

def getGridPartType(gridpart, **parameters):
    """Return the gridpart type (using a function from database.py).
    """
    return myGenerator.getTypeName(gridpart, **parameters)

def get(gridpart, dfmodule, **parameters):
    """Create a gridpart module using the gridpart-database.

    This function creates a python module called gridpart_xxx.py where xxx is a
    number derived from the gridpart type.
    It does this by fetching the typedef and includes from the gridpart-database
    using the given arguments.

    Returns:
        module: the newly created gridpart module

    """
    module = myGenerator.getModule(gridpart, gf=dfmodule._typeName, extra_includes=dfmodule._includes, **parameters)
    setattr(module.GridPart, "_module", module)
    setattr(module.GridPart, "interpolate", interpolate )
    return module

def create(gridpart, gf, **parameters):
    """Get a GridPart.

    Call get() and create a C++ gridpart class

    Args:
        gridpart (string): the identifier for the scheme type to use
        gf; main argument, will in some form conatin the host grid part
    Returns:
        GridPart: the constructed GridPart
    """
    module = get(gridpart, gf._module, **parameters)
    return module.GridPart(gf)

#############################################
if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
