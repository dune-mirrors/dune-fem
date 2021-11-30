from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(space, vec, name="tmp", **unused):
    """create a discrete function - using the fem numpy storage as linear algebra backend

    Args:
        space: discrete space

    Returns:
        DiscreteFunction: the constructed discrete function
    """

    includes = [ "dune/fem/function/vectorfunction/managedvectorfunction.hh", "dune/python/common/numpyvector.hh" ] + space._module.cppIncludes
    spaceType = space._module.cppTypeName
    field = space.field
    typeName = "Dune::Fem::VectorDiscreteFunction< " +\
          spaceType + ", Dune::FemPy::NumPyVector< " + field + " > >"

    return module("numpy", includes, typeName).DiscreteFunction(space,name,vec)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
