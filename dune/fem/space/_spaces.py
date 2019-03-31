from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
from dune.generator import Constructor, Method
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

# dimrange parameter in space creation is deprecated!
def checkDeprecated_dimrange( dimRange, dimrange ):
    import warnings
    if dimRange is not None:
        return dimRange
    elif dimrange is not None:
        warnings.warn('DiscreteFunctionSpace: parameter dimrange for spaces is deprecated, use dimRange instead!\n')
        return dimrange
    else:
        return None

def dgonb(gridView, order=1, dimRange=None, field="double", storage=None,
          interiorQuadratureOrders=None, skeletonQuadratureOrders=None,
          scalar=False, dimrange=None):
    """create a discontinous galerkin space with elementwise orthonormal basis functions

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage, codegen

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")

    if dimRange < 1:
        raise KeyError(\
            "Parameter error in DiscontinuosGalerkinSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in DiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = gridView._includes + [ "dune/fem/space/discontinuousgalerkin.hh" ]
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::DiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    if interiorQuadratureOrders is not None or skeletonQuadratureOrders is not None:
        codegen(spc,interiorQuadratureOrders,skeletonQuadratureOrders)
        typeName = "Dune::Fem::DiscontinuousGalerkinSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
          "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + ", " +\
          "Dune::Fem::CodegenStorage" +\
          " >"
        spc = module(field, includes, typeName,
                    interiorQuadratureOrders=interiorQuadratureOrders,
                    skeletonQuadratureOrders=skeletonQuadratureOrders,storage=storage,
                    scalar=scalar,
                    ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()

def dgonbhp(gridView, order=1, dimRange=None, field="double", storage=None, scalar=False, dimrange=None):
    """create a discontinous galerkin space with elementwise orthonormal basis functions capable of hp-adaptation

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in hpDG::OrthogonalDiscontinuosGalerkinSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in hpDG::OrthogonalDiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/hpdg/orthogonal.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()

def dglegendre(gridView, order=1, dimRange=None, field="double", storage=None, hierarchical=True, scalar=False, dimrange=None):
    """create a discontinous galerkin space with elementwise legendre tensor product basis function

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend
        hierarchical: is true sort the shape function according to their
                      polynomial order otherwise a product arrangement is used

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")

    # if not (len(gridView.indexSet.types(0)) == 1 and
    #         gridView.indexSet.types(0)[0].isCube):
    if not (gridView.type.isCube):
        raise KeyError(\
            "the `dglegendre' space can only be used with a fully "+
            "quadrilateral grid")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in DiscontinuosGalerkinSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in DiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    className = "Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace" \
                if hierarchical else \
                "Dune::Fem::LegendreDiscontinuousGalerkinSpace"
    typeName = className + "< "\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()

def dglegendrehp(gridView, order=1, dimRange=None, field="double", storage=None, scalar=False, dimrange=None):
    """create a discontinous galerkin space with elementwise legendre tensor product basis function capable of hp-adaptation

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    # if not (len(gridView.indexSet.types(0)) == 1 and
    #         gridView.indexSet.types(0)[0].isCube):
    if not (gridView.type.isCube):
        raise KeyError(\
            "the `dglegendrehp' space can only be used with a fully "+
            "quadrilateral grid")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in DiscontinuosGalerkinSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in hpDG::HierarchicLegendreDiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/hpdg/legendre.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()

def dglagrange(gridView, order=1, dimRange=None, field="double", storage=None,
            interiorQuadratureOrders=None, skeletonQuadratureOrders=None,
            scalar=False, dimrange=None):
    """create a discontinous galerkin space with elementwise lagrange basis function

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage, codegen

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in LagrangeDiscontinuosGalerkinSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in LagrangeDiscontinuousGalerkinSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    dimw = gridView.dimWorld

    includes = gridView._includes + [ "dune/fem/space/discontinuousgalerkin.hh" ]
    typeName = "Dune::Fem::DGLagrangeSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " > >"
    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView,order])
    if interiorQuadratureOrders is not None or skeletonQuadratureOrders is not None:
        codegen(spc,interiorQuadratureOrders,skeletonQuadratureOrders)
        typeName = "Dune::Fem::LagrangeDiscontinuousGalerkinSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
          "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + ", " +\
          "Dune::Fem::CodegenStorage" +\
          " >"
        spc = module(field, includes, typeName,
                    interiorQuadratureOrders=interiorQuadratureOrders,
                    skeletonQuadratureOrders=skeletonQuadratureOrders,
                    storage=storage,
                    scalar=scalar,
                    ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()


def lagrange(gridView, order=1, dimRange=None, field="double", storage=None,
            interiorQuadratureOrders=None, skeletonQuadratureOrders=None,
            scalar=False, dimrange=None):
    """create a Lagrange space

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage, codegen

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 1")
    if field == "complex":
        field = "std::complex<double>"

    includes = gridView._includes + [ "dune/fem/space/lagrange.hh" ]
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::LagrangeSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " > >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView,order])
    if interiorQuadratureOrders is not None or skeletonQuadratureOrders is not None:
        codegen(spc,interiorQuadratureOrders,skeletonQuadratureOrders)
        typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
          "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + ", " +\
          "Dune::Fem::CodegenStorage" +\
          " >"
        spc = module(field, includes, typeName,
                    interiorQuadratureOrders=interiorQuadratureOrders,
                    skeletonQuadratureOrders=skeletonQuadratureOrders,storage=storage,
                    scalar=scalar,
                    ctorArgs=[gridView,order])
    return spc.as_ufl()

def lagrangehp(gridView, maxOrder=1, dimRange=None, field="double", storage=None,
               interiorQuadratureOrders=None, skeletonQuadratureOrders=None,
               scalar=False, dimrange=None):
    """create a Lagrange space

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage, codegen

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if maxOrder < 1:
        raise KeyError(\
            "Parameter error in LagrangeHP with "+
            "maxOrder=" + str(maxOrder) + ": " +\
            "maximum order has to be greater or equal to 1")
    if field == "complex":
        field = "std::complex<double>"

    includes = gridView._includes + [ "dune/fem/space/padaptivespace/lagrange.hh" ]
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::PAdaptiveLagrangeSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(maxOrder) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    if interiorQuadratureOrders is not None or skeletonQuadratureOrders is not None:
        codegen(spc,interiorQuadratureOrders,skeletonQuadratureOrders)
        typeName = "Dune::Fem::PAdaptiveLagrangeSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
          "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(maxOrder) + ", " +\
          "Dune::Fem::CodegenStorage" +\
          " >"
        spc = module(field, includes, typeName,
                    interiorQuadratureOrders=interiorQuadratureOrders,
                    skeletonQuadratureOrders=skeletonQuadratureOrders,storage=storage,
                    scalar=scalar,
                    ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()

def finiteVolume(gridView, dimRange=None, field="double", storage=None, scalar=False, dimrange=None):
    """create a finite volume space

    A finite volume space is a discontinuous function space, using the element
    indicator function as its sole basis function.

    Args:
        gridView: the underlying grid view
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    #if dimrange is not None:
    #    raise ValueError('parameter dimrange for spaces is deprecated, use dimRange instead!')

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError("invalid dirRange: " + str(dimRange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = ["dune/fem/space/finitevolume.hh" ] + gridView._includes
    functionSpaceType = "Dune::Fem::FunctionSpace< double, " + field + ", " + str(gridView.dimWorld) + ", " + str(dimRange) + " >"
    typeName = "Dune::Fem::FiniteVolumeSpace< " + functionSpaceType + ", Dune::FemPy::GridPart< " + gridView._typeName + " > >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()


def p1Bubble(gridView, dimRange=None, field="double", order=1, storage=None, scalar=False, dimrange=None):
    """create a P1 space enriched with element bubble functions

    Args:
        gridView: the underlying grid part
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in P1BubbleSpace with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if not order == 1:
        raise KeyError(\
            "Parameter error in P1BubbleSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/p1bubble.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::BubbleElementSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " > >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()


def combined(*spaces, **kwargs):
    """create a discrete function space from a tuple of discrete function spaces

    Args:
        spaces: tuple of discrete function spaces

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module, addStorage

    scalar = kwargs.get("scalar",False)

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

    spc = module(combinedField, includes, typeName, constructor,
            storage=combinedStorage,
            scalar=scalar,
            ctorArgs=[spaces])
    try:
        spc.componentNames = kwargs["components"]
    except KeyError:
        pass
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

    scalar = kwargs.get("scalar",False)

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

    storage = lambda _: [None,combinedIncludes+["dune/fem/function/tuplediscretefunction.hh"],
                        "Dune::Fem::TupleDiscreteFunction< " + ", ".join(s.storage[2] for s in spaces) + " >",
                        None,None,None]
    spc = module(combinedField, includes, typeName, constructor, storage=storage,
            scalar=scalar,
            ctorArgs=[spaces])
    try:
        spc.componentNames = kwargs["components"]
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
        print(name,space.componentNames)
        # try:
        df = tupleDiscreteFunction(space, name=name, components=space.componentNames,**kwargs)
        # except AttributeError:
        #     df = tupleDiscreteFunction(space, name=name, **kwargs)
        df.interpolate(func)
        return df
    setattr(spc, "interpolate", lambda *args,**kwargs: interpolate(spc,*args,**kwargs))

    # there is no obvious operator associated with the TupleDF used for this space
    addStorage(spc, lambda _: [None,combinedIncludes+["dune/fem/function/tuplediscretefunction.hh"],
                        "Dune::Fem::TupleDiscreteFunction< " + ", ".join(s.storage[2] for s in spaces) + " >",
                        None,None,None] )
    return spc.as_ufl()

def bdm(gridView, order=1, field="double", storage=None, scalar=False, dimrange=None):
    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if order < 1:
        raise KeyError(\
            "Parameter error in BDMSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1 or 2")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/brezzidouglasmarini.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::BrezziDouglasMariniSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimw) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    return spc.as_ufl()

def raviartThomas(gridView, order=1, field="double", storage=None, scalar=False, dimrange=None):
    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if order > 2:
        raise KeyError(\
            "Parameter error in RTSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 0,1 or 2")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/raviartthomas.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::RaviartThomasSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimw) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " >, " + str(order) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    return spc.as_ufl()

def rannacherTurek(gridView, dimRange=None, field="double", storage=None, scalar=False, dimrange=None):
    from dune.fem.space import module, addStorage

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError(\
                "trying to set up a scalar space with dimRange = " +\
                str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError("invalid dimRange: " + str(dimRange) + " (must be >= 1)")
    if field == "complex":
        field = "std::complex<double>"

    includes = [ "dune/fem/space/second.hh" ] + gridView._includes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::RannacherTurekDiscreteFunctionSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView._typeName + " > >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar,
            ctorArgs=[gridView])
    # addStorage(spc, storage)
    return spc.as_ufl()
