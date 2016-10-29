from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

def dgonb(gridview, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a discontinous galerkin space with elementwise orthonormal basis functions

    Args:
        gridview: the underlying grid view
        order: polynomial order of the finite element functions
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module
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

    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridview._includes
    dimw = gridview.dimWorld
    typeName = "Dune::Fem::DiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(gridview)

def lagrange(view, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a Lagrange space

    Args:
        view: the underlying grid view
        order: polynomial order of the finite element functions
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module
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

    includes = [ "dune/fem/space/lagrange.hh" ] + view._includes
    dimw = view.dimWorld
    typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + view._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(view)

def p1Bubble(gridview, dimrange=1, field="double", order=1, storage=None, **unused):
    """create a P1 space enriched with element bubble functions

    Args:
        gridview: the underlying grid part
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module
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

    includes = [ "dune/fem/space/p1bubble.hh" ] + gridview._includes
    dimw = gridview.dimWorld
    typeName = "Dune::Fem::BubbleElementSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " > >"

    return module(field, storage, includes, typeName).Space(gridview)
