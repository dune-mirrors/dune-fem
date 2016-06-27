from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

from ..generator import generator

myGenerator = generator.Generator("DiscreteFunction")

def getDiscreteFunctionType(df, **parameters):
    """Return the discrete function type (using a function from database.py).
    """
    return myGenerator.getTypeName(df, **parameters)

def get(df, spaceModule, **parameters):
    """Create a discrete function module using the discretefunction-database.

    This function creates a python module called discretefunction_xxx.py where xxx is a
    number derived from the discrete function type.
    It does this by fetching the typedef and includes from the
    discrerefunction-database
    using the given arguments.

    Args:
        df (string): the identifier for the discrete function type to use
        spaceModule (module): space module to build the discrete function on
        parameters (kwargs): parameters used for fixing the discrete function type

            * dimrange=int: dimension of the range type
            * polorder=int: polynomial order of the discrete function

    Returns:
        module: the newly created discrete function module

    """
    module = myGenerator.getModule(df, extra_includes=spaceModule._includes, space=spaceModule._typeName, **parameters)
    setattr(module.DiscreteFunction, "_module", module)
    setattr(module.DiscreteFunction, "_storage", module._selector.grid)
    return module

def create(df, space, **parameters):
    """Get a discrete function

    Call get() and create a C++ discrete function class.

    Notes:
        This is equivalent to::

            discreteFunctionModule = get(df,space._module,parameters)
            df = discreteFunctionModule.DiscreteFunction(space)

    Args:
        df (string): the identifier for the discrete function type to use
        space (Space): space to build the discrete function on
        parameters (kwargs): parameters used for fixing the space type

            * name=str: name of discrete function for output

    Returns:
        DiscreteFunction: the constructed discrete function

    """
    module=get(df, space._module, **parameters)
    name = parameters.get("name","default")
    return module.DiscreteFunction(space,name)


#############################################

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
