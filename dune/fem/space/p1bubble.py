from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(gridpart, dimrange=1, field="double", order=1, storage="adaptive", **unused):
    """create a P1 space enriched with element bubble functions

    Args:
        gridpart: the underlying grid part
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    if dimrange < 1:
        raise KeyError(\
            "Parameter error in P1BubbleSpace with "+
            "dimrange=" + str(dimrange) + ": " +\
            "dimrange has to be greater or equal to 1")
    if not order == 1:
        raise KeyError(\
            "Parameter error in P1BubbleSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/p1bubble.hh" ] + gridpart._module._includes
    gridPart = gridpart._module._typeName
    dimw = gridpart.dimWorld
    typeName = "Dune::Fem::BubbleElementSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      gridPart + " >"

    return module(field, includes, typeName).Space(gridpart)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
