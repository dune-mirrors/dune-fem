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

def _checkDimRangeScalarOrderField(dimRange, scalar, order, field, methodName = None):
    """ method checking parameters for all spaces """

    if methodName is None:
        methodName = sys._getframe().f_back.f_code.co_name

    if dimRange is None:
        dimRange = 1
        scalar = True

    if dimRange > 1 and scalar:
        raise KeyError("In '" + methodName + "': trying to set up a scalar space with dimRange = " +str(dimRange) + ">1")
    if dimRange < 1:
        raise KeyError(\
            "Parameter error in '" + methodName + "' with "+
            "dimRange=" + str(dimRange) + ": " +\
            "dimRange has to be greater or equal to 1")
    if order < 0:
        raise KeyError(\
            "Parameter error in '" + methodName + "' with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 0")
    if field == "complex":
        field = "std::complex<double>"

    return dimRange,scalar,field

# maxOrder parameter in space creation is deprecated!
def checkDeprecated_maxOrder( order, maxOrder ):
    import warnings
    if maxOrder is not None:
        warnings.warn('DiscreteFunctionSpace: parameter maxOrder for spaces is deprecated, use order instead!\n')
        return maxOrder
    else:
        return order

def storageType(codegen):
    if codegen is True or codegen == "codegen":
        return "Dune::Fem::CodegenStorage"
    elif codegen is False or codegen == "caching":
        return "Dune::Fem::CachingStorage"
    elif codegen == "simple":
        return "Dune::Fem::SimpleStorage"
    else:
        raise KeyError(\
            "Parameter error in space construction with "+
            "codege=" + codegen + ": " +\
            "codegen needs to be in [True,False,'codegen','caching','simple'")

def dgonb(gridView, order=1, dimRange=None, field="double",
          storage=None, scalar=False, dimrange=None, codegen=True):
    """create a discontinuous galerkin space with elementwise orthonormal basis functions

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend
        caching: true if basis functions are cached for each quadrature

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = gridView.cppIncludes + [ "dune/fem/space/discontinuousgalerkin.hh" ]
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::DiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + ", " + storageType(codegen) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])

    return spc.as_ufl()

def dgonbhp(gridView, order=1, dimRange=None, field="double",
            storage=None, scalar=False, dimrange=None, codegen=True):
    """create a discontinuous galerkin space with elementwise orthonormal basis functions capable of hp-adaptation

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = [ "dune/fem/space/hpdg/orthogonal.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::hpDG::OrthogonalDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + " ," + storageType(codegen) + " >"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()

def dglegendre(gridView, order=1, dimRange=None, field="double",
               storage=None, hierarchical=True, scalar=False, dimrange=None, codegen=True):
    """create a discontinuous galerkin space with elementwise legendre tensor product basis function

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

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    if not (gridView.type.isCube):
        raise KeyError(\
            "the `dglegendre' space can only be used with a fully "+
            "quadrilateral/hexahedral grid")

    includes = [ "dune/fem/space/discontinuousgalerkin.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    className = "Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace" \
                if hierarchical else \
                "Dune::Fem::LegendreDiscontinuousGalerkinSpace"
    typeName = className + "< "\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + ", "+\
      storageType(codegen) + ">"
    ctorArgs = [gridView]
    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=ctorArgs)
    return spc.as_ufl()

###########################################################################
##  Legendre hpDG space
###########################################################################
def dglegendrehp(gridView, order=1, dimRange=None, field="double",
                 storage=None, scalar=False, dimrange=None, codegen=True):
    """create a discontinuous galerkin space with elementwise legendre tensor product basis function capable of hp-adaptation

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    if not (gridView.type.isCube):
        raise KeyError(\
            "the `dglegendrehp' space can only be used with a fully "+
            "quadrilateral grid")

    includes = [ "dune/fem/space/hpdg/legendre.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::hpDG::HierarchicLegendreDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + ", " +\
      storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()

###########################################################################
##  ansiotropic hpdg space
###########################################################################
def dganisotropic(gridView, order=1, dimRange=None, field="double",
                  storage=None, scalar=False, dimrange=None, codegen=True):
    """create a discontinuous Galerkin space with element wise anisotropic (Legendre)
       basis function capable of hp-adaptation

    Args:
        gridView: the underlying grid view
        order: maximal polynomial order or list of orders for each spatial
               direction of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # obtain max order
    maxOrder = order
    isAniso = False
    if not isinstance(order,int) and (isinstance(order,list) or isinstance(order,tuple)):
        assert len(order) == gridView.dimGrid, "If order is provided as list or tuple the length needs to match the dimension of the grid!"
        maxOrder = max(order)
        isAniso = True

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, maxOrder, field)

    includes = [ "dune/fem/space/hpdg/anisotropic.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::hpDG::AnisotropicDiscontinuousGalerkinSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(maxOrder) + ", " + storageType(codegen) + ">"

    ## constructor taking vector with orders
    constructor = Constructor(['pybind11::object gridView', 'std::vector<int> orders'],
                              ['return new DuneType( Dune::FemPy::gridPart< typename DuneType::GridPartType::GridViewType >( gridView ), orders );'],
                              ['"gridView"_a', '"orders"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()'])

    spc = module(field, includes, typeName, constructor, storage=storage,
                 scalar=scalar, codegen=codegen,
                 ctorArgs = [gridView, order] if isAniso else [gridView])
    return spc.as_ufl()

def dglagrange(gridView, order=1, dimRange=None, field="double", storage=None,
               scalar=False, dimrange=None, codegen=True, pointType=None):
    """create a discontinuous galerkin space with elementwise lagrange basis function

    Args:
        gridView: the underlying grid view
        order: polynomial order of the finite element functions
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    dimw = gridView.dimWorld

    includes = gridView.cppIncludes + [ "dune/fem/space/discontinuousgalerkin.hh" ]
    if pointType is None:
        typeName = "Dune::Fem::LagrangeDiscontinuousGalerkinSpace< " +\
          "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
          "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + " ," +\
          storageType(codegen) + ">"
        ctorArgs=[gridView]
    else:
        def errorText(key):
            return f"the `dglagrange({key}) space can only be used with a fully "+ \
                    "quadrilateral/hexahedral grid"

        includes += ["dune/fem/space/localfiniteelement/quadratureinterpolation.hh"]
        if pointType.lower() == "equidistant":
            pointSet = 'Dune::EquidistantPointSetDerived'
        elif pointType.lower() == "lobatto":
            if not (gridView.type.isCube):
                raise KeyError(errorText("lobatto"))
            pointSet = 'Dune::GaussLobattoPointSet'
        elif pointType.lower() == "gauss":
            if not (gridView.type.isCube):
                raise KeyError(errorText("gauss"))
            pointSet = 'Dune::GaussLegendrePointSet'
        elif pointType.lower() == "cellcenters":
            if not (gridView.type.isCube):
                raise KeyError(errorText("cellcenters"))
            pointSet = 'Dune::CellCentersPointSet'
        else:
            raise KeyError(
                "Parameter error in LagrangeDiscontinuousGalerkinSpace with point set type " +
                pointType + "not known.")
        if False: # none fixed order version
            typeName = "Dune::Fem::DGLagrangeSpace< " +\
              "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
              "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, "+\
                  pointSet +\
              ", Dune::Fem::CachingStorage >"
            ctorArgs=[gridView,order]
        else: # fixed order version
            typeName = "Dune::Fem::FixedOrderDGLagrangeSpace< " +\
              "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
              "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >,"+\
              str(order) + ", " + pointSet + ", " +\
              storageType(codegen) + ">"
            ctorArgs=[gridView]

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=ctorArgs)
    return spc.as_ufl()

def dglagrangelobatto(*args, **kwargs):
    """create a discontinuous galerkin space with elementwise lagrange basis function
       and Legendre-Gauss-Lobatto interpolation points.

    Args:
        same as dglagrange

    Returns:
        Space: the constructed Space, same as dglagrange(pointType='lobatto')
    """
    return dglagrange(*args, **kwargs, pointType='lobatto')

def lagrange(gridView, order=1, dimRange=None, field="double", storage=None,
             scalar=False, dimrange=None, codegen=True):
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

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if order < 1:
        raise KeyError(\
            "Parameter error in LagrangeSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be greater or equal to 1")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = gridView.cppIncludes + [ "dune/fem/space/lagrange.hh" ]
    dimw = gridView.dimWorld
    # for order equal or lesser than 6 we can use
    # DynamicLagrangeDiscreteFunctionSpace to avoid re-compilation
    if order <= 6:
        typeName = "Dune::Fem::DynamicLagrangeDiscreteFunctionSpace< " +\
                   "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
                   "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >"
    else:
        typeName = "Dune::Fem::LagrangeDiscreteFunctionSpace< " +\
                   "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
                   "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order)

    spcTypeName = typeName + ", " + storageType(codegen) + ">"
    spc = module(field, includes, spcTypeName, storage=storage,
                 scalar=scalar, codegen=codegen,
                 ctorArgs=[gridView, order])

    return spc.as_ufl()

def lagrangehp(gridView, order=1, dimRange=None, field="double", storage=None,
               scalar=False, maxOrder=None, dimrange=None, codegen=True):
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

    from dune.fem.space import module

    order    = checkDeprecated_maxOrder( order=order, maxOrder=maxOrder )
    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if order < 1:
        raise KeyError(\
            "Parameter error in LagrangeHP with "+
            "order=" + str(order) + ": " +\
            "maximum order has to be greater or equal to 1")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    # set maxOrder of space to 6 even though used maxOrder given by order
    # may be smaller, this way compilation time can be reduced
    maxOrder = 6 if order <= 6 else order

    includes = gridView.cppIncludes + [ "dune/fem/space/padaptivespace/lagrange.hh" ]
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::PAdaptiveLagrangeSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(maxOrder) + ", " +\
      storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView,order])
    return spc.as_ufl()

def finiteVolume(gridView, dimRange=None, field="double",
                 storage=None, scalar=False, dimrange=None,
                 codegen='simple'):
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

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, 0, field)

    includes = ["dune/fem/space/finitevolume.hh" ] + gridView.cppIncludes
    functionSpaceType = "Dune::Fem::FunctionSpace< double, " + field + ", " + str(gridView.dimWorld) + ", " + str(dimRange) + " >"
    typeName = "Dune::Fem::FiniteVolumeSpace< " + functionSpaceType +\
               ", Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, 0 , " +\
               storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()


def p1Bubble(gridView, dimRange=None, field="double", order=1,
             storage=None, scalar=False, dimrange=None, codegen=True):
    """create a P1 space enriched with element bubble functions

    Args:
        gridView: the underlying grid part
        dimRange: dimension of the range space
        field: field of the range space
        storage: underlying linear algebra backend

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if not order == 1:
        raise KeyError(\
            "Parameter error in P1BubbleSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = [ "dune/fem/space/p1bubble.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::BubbleElementSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " > >, " +\
      storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()


def combined(*spaces, **kwargs):
    """create a discrete function space from a tuple of discrete function spaces

    Args:
        spaces: tuple of discrete function spaces

    Returns:
        Space: the constructed Space
    """

    from dune.fem.space import module

    scalar = kwargs.get("scalar",False)
    codegen = kwargs.get("codegen",True)

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
        includes += space.cppIncludes
    typeName = "Dune::Fem::TupleDiscreteFunctionSpace< " + ", ".join([space.cppTypeName for space in spaces]) + " >"

    constructor = Constructor(['typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple'],
                              ['return new DuneType( spaceTuple);'],
                              ['"spaceTuple"_a', 'pybind11::keep_alive<1,2>()'])

    spc = module(combinedField, includes, typeName, constructor,
            storage=combinedStorage,
            scalar=scalar, codegen=codegen,
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

    scalar  = kwargs.get("scalar",False)
    codegen = kwargs.get("codegen",False)

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
        includes += space.cppIncludes
    typeName = "Dune::Fem::TupleDiscreteFunctionSpace< " + ", ".join([space.cppTypeName for space in spaces]) + " >"

    constructor = Constructor(['typename DuneType::DiscreteFunctionSpaceTupleType spaceTuple'],
                              ['return new DuneType( spaceTuple);'],
                              ['"spaceTuple"_a', 'pybind11::keep_alive<1,2>()'])

    storage = lambda _: [None,combinedIncludes+["dune/fem/function/tuplediscretefunction.hh"],
                        "Dune::Fem::TupleDiscreteFunction< " + ", ".join(s.storage[2] for s in spaces) + " >",
                        None,None,None]
    spc = module(combinedField, includes, typeName, constructor, storage=storage,
            scalar=scalar, codegen=codegen,
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

def bdm(gridView, order=1, dimRange=None,
        field="double", storage=None, scalar=False, dimrange=None, codegen=True):
    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if order < 1:
        raise KeyError(\
            "Parameter error in BDMSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 1 or 2")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = [ "dune/fem/space/brezzidouglasmarini.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::BrezziDouglasMariniSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimw) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + " ," +\
      storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()

def raviartThomas(gridView, order=1, dimRange=None,
                  field="double", storage=None, scalar=False, dimrange=None, codegen=True):
    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if order > 2:
        raise KeyError(\
            "Parameter error in RTSpace with "+
            "order=" + str(order) + ": " +\
            "order has to be equal to 0,1 or 2")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, order, field)

    includes = [ "dune/fem/space/raviartthomas.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::RaviartThomasSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimw) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " >, " + str(order) + " ," +\
      storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()

def rannacherTurek(gridView, dimRange=None,
                   field="double", storage=None, scalar=False, dimrange=None, codegen=True):
    from dune.fem.space import module

    dimRange = checkDeprecated_dimrange( dimRange=dimRange, dimrange=dimrange )

    if dimRange and dimRange < 1:
        raise KeyError("invalid dimRange: " + str(dimRange) + " (must be >= 1)")

    # check requirements on parameters
    dimRange, scalar, field = _checkDimRangeScalarOrderField(dimRange, scalar, 1, field)

    includes = [ "dune/fem/space/rannacherturek.hh" ] + gridView.cppIncludes
    dimw = gridView.dimWorld
    typeName = "Dune::Fem::RannacherTurekSpace< " +\
      "Dune::Fem::FunctionSpace< double, " + field + ", " + str(dimw) + ", " + str(dimRange) + " >, " +\
      "Dune::FemPy::GridPart< " + gridView.cppTypeName + " > ," + storageType(codegen) + ">"

    spc = module(field, includes, typeName, storage=storage,
            scalar=scalar, codegen=codegen,
            ctorArgs=[gridView])
    return spc.as_ufl()
