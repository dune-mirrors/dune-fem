from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(gridpart, order=1, dimrange=1, field="double", storage="adaptive", **unused):
    """create a discontinous galerkin space with elementwise orthonormal basis functions

    Args:
        gridpart: the underlying grid part
        order: polynomial order of the finite element functions
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    if dimrange < 1:
        raise KeyError(\
            "Parameter error in DiscontinuosGalerkinSpace with "+
            "dimrange=" + str(dimrange) + ": " +\
            "dimrange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in DiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridpart._module._includes
    gridPart = gridpart._module._typeName
    dimw = gridpart.dimWorld
    typeName = "Dune::Fem::DiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      gridPart + ", " + str(order) + " >"

    return module(field, includes, typeName).Space(gridpart)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
