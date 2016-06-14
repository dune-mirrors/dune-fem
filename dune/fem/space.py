from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

from ..generator import generator

myGenerator = generator.Generator("Space")

def getSpaceType(space, **parameters):
    """Return the space type (using a function from database.py).
    """
    return myGenerator.getTypeName(space, **parameters)

def get(space, gridModule, **parameters):
    """Create a space module using the space-database.

    This function creates a python module called space_xxx.py where xxx is a
    number derived from the space type.
    It does this by fetching the typedef and includes from the space-database
    using the given arguments.

    Args:
        space (string): the identifier for the space type to use
        gridModule (module): grid module to build the space on
        parameters (kwargs): parameters used for fixing the space type

            * dimrange=int: dimension of the range type
            * polorder=int: polynomial order of the space

    Returns:
        module: the newly created space module

    """
    return myGenerator.getModule(space, extra_includes=gridModule._includes, gridpart=gridModule._typeName, **parameters)

def create(space, grid, **parameters):
    """Get a Space

    Call get() and create a C++ space class.

    Notes:
        This is equivalent to::

            spaceModule = get(space,grid._module,parameters)
            space = spaceModule.Space(grid)

    Args:
        space (string): the identifier for the space type to use
        grid (LeafGrid): grid to build the space on
        parameters (kwargs): parameters used for fixing the grid type

            * dimrange=int: dimension of the range type
            * polorder=int: polynomial order of the space

    Returns:
        Space: the constructed space

    """
    module=get(space, grid._module, **parameters)
    return module.Space(grid)


#############################################

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
