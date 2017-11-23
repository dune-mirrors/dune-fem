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

def dgonbhp(gridview, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a discontinous galerkin space with elementwise orthonormal basis functions capable of hp-adaptation

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
            "Parameter error in hpDG::OrthogonalDiscontinuosGalerkinSpace with "+
            "dimrange=" + str(dimrange) + ": " +\
            "dimrange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in hpDG::OrthogonalDiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/hpdg/orthogonal.hh" ] + gridview._includes
    dimw = gridview.dimWorld
    typeName = "Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(gridview)

def dglegendre(gridview, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a discontinous galerkin space with elementwise legendre tensor product basis function

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
    typeName = "Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(gridview)

def dglegendrehp(gridview, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a discontinous galerkin space with elementwise legendre tensor product basis function capable of hp-adaptation

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
            "Parameter error in hpDG::HierarchicLegendreDiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/hpdg/legendre.hh" ] + gridview._includes
    dimw = gridview.dimWorld
    typeName = "Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(gridview)

def dglagrange(gridview, order=1, dimrange=1, field="double", storage=None, **unused):
    """create a discontinous galerkin space with elementwise lagrange basis function

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

    dimw = gridview.dimWorld

    # includes = [ "dune/fem/space/lagrange.hh" ] + gridview._includes
    # typeName = "Dune::Fem::DGLagrangeSpace< " +\
    #   "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
    #   "Dune::FemPy::GridPart< " + gridview._typeName + " > >"
    # return module(field, storage, includes, typeName).Space(gridview,order)
    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridview._includes
    typeName = "Dune::Fem::LagrangeDiscontinuousGalerkinSpace< " +\
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


def finiteVolume(gridview, dimrange=1, field="double", storage=None, **unused):
    """create a finite volume space

    A finite volume space is a discontinuous function space, using the element
    indicator function as its sole basis function.

    Args:
        gridview: the underlying grid view
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module
    if dimrange < 1:
        raise KeyError("invalid dimRange: " + str(dimrange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = ["dune/fem/space/finitevolume.hh" ] + gridview._includes
    functionSpaceType = "Dune::Fem::FunctionSpace< double, " + field + ", " + str(gridview.dimWorld) + ", " + str(dimrange) + " >"
    typeName = "Dune::Fem::FiniteVolumeSpace< " + functionSpaceType + ", Dune::FemPy::GridPart< " + gridview._typeName + " > >"

    return module(field, storage, includes, typeName).Space(gridview)


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


def combined(*spaces, **unused):
    """create a discrete function space from a tuple of discrete function spaces

    Args:
        spaces: tuple of discrete function spaces

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    if not spaces:
        raise Exception("Cannot create TupleDiscreteFunctionSpace from empty tuple of discrete function spaces")
    combinedStorage = None
    combinedField = None
    for space in spaces:
        storage, _, _, _, _ = space.storage
        if combinedStorage and (combinedStorage != storage):
            raise Exception("Cannot create TupleDiscreteFunctionSpace with different types of storage")
        else:
            combinedStorage = storage
        if combinedField and (combinedField != space.field):
            raise Exception("Cannot create TupleDiscreteFunctionSpace with different field types")
        else:
            combinedField = space.field

    includes = ["dune/fem/space/combinedspace/tuplespace.hh"]
    for space in spaces:
        includes += space._includes
    typeName = "Dune::Fem::TupleDiscreteFunctionSpace< " + ", ".join([space._typeName for space in spaces]) + " >"

    return module(combinedField, combinedStorage, includes, typeName).Space(spaces[0].grid)

def bdm(view, order=1, field="double", storage=None, **unused):
    from dune.fem.space import module
    if order < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1 or 2")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/brezzidouglasmarini.hh" ] + view._includes
    dimw = view.dimWorld
    typeName = "Dune::Fem::BDMDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimw) + " >, " +\
      "Dune::FemPy::GridPart< " + view._typeName + " >, " + str(order) + " >"

    return module(field, storage, includes, typeName).Space(view)

def rannacherTurek(view, dimrange=1, field="double", storage=None, **unused):
    from dune.fem.space import module
    if dimrange < 1:
        raise KeyError("invalid dimRange: " + str(dimrange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/second.hh" ] + view._includes
    dimw = view.dimWorld
    typeName = "Dune::Fem::RannacherTurekDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + view._typeName + " > >"

    return module(field, storage, includes, typeName).Space(view)
