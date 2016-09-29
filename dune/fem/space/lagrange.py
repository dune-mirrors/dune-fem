from __future__ import absolute_import, division, print_function, unicode_literals

from . import module

def create(gridview, order=1, dimrange=1, field="double", storage="adaptive", **unused):
    """create a Lagrange space

    Args:
        gridview: the underlying grid view
        order: polynomial order of the finite element functions
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    if dimrange < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "dimrange=" + str(dimrange) + ": " +\
            "dimrange has to be greater or equal to 1")
    if order < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 1")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/lagrange.hh" ] + gridview._module._includes
    dimw = gridview.dimWorld
    typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._module._typeName + " >, " + str(order) + " >"

    return module(field, includes, typeName).Space(gridview)

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
