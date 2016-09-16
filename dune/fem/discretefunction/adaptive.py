from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(space, name="tmp", **unused):
    """create a discrete function - using the fem adaptive storage as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    includes = [ "dune/fem/function/adaptivefunction.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "Dune::Fem::AdaptiveDiscreteFunction< " + spaceType + " >"

    return module("fem", includes, typeName).DiscreteFunction(space,name)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
