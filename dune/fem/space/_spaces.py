from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
from dune.generator import Constructor, Method
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

    from dune.fem.space import module, addStorage
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

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()

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

    from dune.fem.space import module, addStorage
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

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()

def dglegendre(gridview, order=1, dimrange=1, field="double", storage=None, hierarchical=True, **unused):
    """create a discontinous galerkin space with elementwise legendre tensor product basis function

    Args:
        gridview: the underlying grid view
        order: polynomial order of the finite element functions
        dimrange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend
        hierarchical: is true sort the shape function according to their
                      polynomial order otherwise a product arrangement is used

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage
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
    className = "Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace" \
                if hierarchical else \
                "Dune::Fem::LegendreDiscontinuousGalerkinSpace"
    typeName = className + "< "\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()

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

    from dune.fem.space import module, addStorage
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

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()

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

    from dune.fem.space import module, addStorage
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

    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridview._includes
    typeName = "Dune::Fem::LagrangeDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridview._typeName + " >, " + str(order) + " >"
    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()


def lagrange(view, order=1, dimrange=1, field="double", storage=None,
            interiorQuadratureOrders=None, skeletonQuadratureOrders=None,
            **unused):
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

    from dune.fem.space import module, addStorage, codegen
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

    includes = view._includes + [ "dune/fem/space/lagrange.hh" ]
    dimw = view.dimWorld
    typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + view._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName).Space(view)
    if interiorQuadratureOrders is not None or skeletonQuadratureOrders is not None:
        codegen(spc,interiorQuadratureOrders,skeletonQuadratureOrders)
        typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
          "Dune::FemPy::GridPart< " + view._typeName + " >, " + str(order) + ", " +\
          "Dune::Fem::CodegenStorage" +\
          " >"
        spc = module(field, includes, typeName,
                    interiorQuadratureOrders=interiorQuadratureOrders,
                    skeletonQuadratureOrders=skeletonQuadratureOrders).Space(view)
    addStorage(spc, storage)
    return spc.as_ufl()

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

    from dune.fem.space import module, addStorage
    if dimrange < 1:
        raise KeyError("invalid dimRange: " + str(dimrange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = ["dune/fem/space/finitevolume.hh" ] + gridview._includes
    functionSpaceType = "Dune::Fem::FunctionSpace< double, " + field + ", " + str(gridview.dimWorld) + ", " + str(dimrange) + " >"
    typeName = "Dune::Fem::FiniteVolumeSpace< " + functionSpaceType + ", Dune::FemPy::GridPart< " + gridview._typeName + " > >"

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()


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

    from dune.fem.space import module, addStorage
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

    spc = module(field, includes, typeName).Space(gridview)
    addStorage(spc, storage)
    return spc.as_ufl()


def combined(*spaces, **kwargs):
    """create a discrete function space from a tuple of discrete function spaces

    Args:
        spaces: tuple of discrete function spaces

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    if not spaces:
        raise Exception("Cannot create TupleDiscreteFunctionSpace from empty tuple of discrete function spaces")
    spaces = tuple([s.__impl__ for s in spaces])
    combinedStorage = None
    combinedField = None
    for space in spaces:
        storage, _, _, _, _, _ = space.storage
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

    constructor = Constructor(['typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple'],
                              ['return new DuneType( spaceTuple);'],
                              ['"spaceTuple"_a', 'pybind11::keep_alive<1,2>()'])

    mod = module(combinedField, includes, typeName, constructor)
    try:
        mod.Space.componentNames = kwargs["components"]
    except KeyError:
        pass
    spc = mod.Space(spaces)
    addStorage(spc, combinedStorage)
    return spc.as_ufl()

def product(*spaces, **kwargs):
    """create a discrete function space from a tuple of discrete function spaces

    Args:
        spaces: tuple of discrete function spaces

    Returns:
        Space: the constructed Space
    """
    from dune.fem.space import module, addStorage
    from dune.fem.function import tupleDiscreteFunction

    if len(spaces)==1 and (isinstance(spaces[0],list) or isinstance(spaces[0],tuple)):
        spaces = spaces[0]
    if not spaces:
        raise Exception("Cannot create TupleDiscreteFunctionSpace from empty tuple of discrete function spaces")

    spaces = tuple([s.__impl__ for s in spaces])

    combinedField = None
    combinedIncludes = None
    for space in spaces:
        storage, _, _, _, _, _ = space.storage
        if combinedField and (combinedField != space.field):
            raise Exception("Cannot create TupleDiscreteFunctionSpace with different field types")
        else:
            combinedField = space.field
    combinedIncludes = [i for s in spaces for i in s.storage[1]]

    includes = ["dune/fem/space/combinedspace/tuplespace.hh"]
    for space in spaces:
        includes += space._includes
    typeName = "Dune::Fem::TupleDiscreteFunctionSpace< " + ", ".join([space._typeName for space in spaces]) + " >"

    constructor = Constructor(['typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple'],
                              ['return new DuneType( spaceTuple);'],
                              ['"spaceTuple"_a', 'pybind11::keep_alive<1,2>()'])
    mod = module(combinedField, includes, typeName, constructor)
    try:
        mod.Space.componentNames = kwargs["components"]
    except KeyError:
        pass
    def interpolate(space, func, name=None, **kwargs):
        """interpolate a function into a discrete function space

        Args:
            space: discrete function space to interpolate into
            func:  function to interpolate
            name:  name of the resulting discrete function

        Returns:
            DiscreteFunction: the constructed discrete function
        """
        if name is None: name = func.name
        try:
            df = tupleDiscreteFunction(space, name=name, components=space.componentNames,**kwargs)
        except AttributeError:
            df = tupleDiscreteFunction(space, name=name, **kwargs)
        df.interpolate(func)
        return df
    setattr(mod.Space, "interpolate", interpolate)

    # there is no obvious operator associated with the TupleDF used for this space
    spc = mod.Space(spaces)
    addStorage(spc, lambda _: [None,combinedIncludes+["dune/fem/function/tuplediscretefunction.hh"],
                        "Dune::Fem::TupleDiscreteFunction< " + ", ".join(s.storage[2] for s in spaces) + " >",
                        None,None,None] )
    return spc.as_ufl()

def bdm(view, order=1, field="double", storage=None, **unused):
    from dune.fem.space import module, addStorage
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

    spc = module(field, includes, typeName).Space(view)
    addStorage(spc, storage)
    return spc.as_ufl()

def rannacherTurek(view, dimrange=1, field="double", storage=None, **unused):
    from dune.fem.space import module, addStorage
    if dimrange < 1:
        raise KeyError("invalid dimRange: " + str(dimrange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/second.hh" ] + view._includes
    dimw = view.dimWorld
    typeName = "Dune::Fem::RannacherTurekDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimrange) + " >, " +\
      "Dune::FemPy::GridPart< " + view._typeName + " > >"

    spc = module(field, includes, typeName).Space(view)
    addStorage(spc, storage)
    return spc.as_ufl()
