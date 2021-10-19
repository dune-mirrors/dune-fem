from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(space, name="tmp", **unused):
    """create a discrete function - using the fem adaptive storage as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    includes = [ "dune/fem/function/adaptivefunction.hh" ] + space._module.cppIncludes
    spaceType = space._module.cppTypeName
    typeName = "Dune::Fem::AdaptiveDiscreteFunction< " + spaceType + " >"

    return module("fem", includes, typeName).DiscreteFunction(space,name)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
