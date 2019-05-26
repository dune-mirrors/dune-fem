from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(space, name="tmp", **unused):
    """create a discrete function - using the eigen library as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    includes = [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/fem/storage/eigenvector.hh" ] + space._module._includes
    spaceType = space._module._typeName
    field = space.field
    typeName = "Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< " +\
            spaceType + ", Dune::Fem::EigenVector< " + field + " > > >"

    return module("eigen", includes, typeName).DiscreteFunction(space,name)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
