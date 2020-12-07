from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
import warnings

from ufl.equation import Equation

from dune.generator import Constructor, Method

logger = logging.getLogger(__name__)

def getSolver(solver,storage,default):
    if not solver:
        return default(storage)
    if isinstance(solver,str):
        return default(storage,solver)
    else:
        tupLen = len(solver)
        import dune.create as create
        if (tupLen > 1):
            return create.solver(solver[0],storage,solver[1])
        else:
            return default(storage,solver[0])

def femscheme(includes, space, solver, operator, modelType):
    storageStr, dfIncludes, dfTypeName, linearOperatorType, defaultSolver, backend = space.storage
    _, solverIncludes, solverTypeName,param = getSolver(solver,space.storage,defaultSolver)

    includes += ["dune/fem/schemes/femscheme.hh"] +\
                space._includes + dfIncludes + solverIncludes +\
                ["dune/fempy/parameter.hh"]
    spaceType = space._typeName
    if modelType is None:
        includes += ["dune/fem/schemes/diffusionmodel.hh"]
        modelType = "DiffusionModel< " +\
              "typename " + spaceType + "::GridPartType, " +\
              spaceType + "::dimRange, " +\
              spaceType + "::dimRange, " +\
              "typename " + spaceType + "::RangeFieldType >"
    operatorType = operator(linearOperatorType, modelType)
    typeName = "FemScheme< " + operatorType + ", " + solverTypeName + " >"
    return includes, typeName

def femschemeModule(space, model, includes, solver, operator, *args,
        parameters={},
        modelType = None, ctorArgs={}):
    from . import module
    _, _, _, _, defaultSolver, backend = space.storage
    _, _, _, param = getSolver(solver,space.storage,defaultSolver)
    includes, typeName = femscheme(includes, space, solver, operator, modelType)
    parameters.update(param)
    mod = module(includes, typeName, *args, backend=backend)
    scheme = mod.Scheme(space, model, parameters=parameters, **ctorArgs)
    scheme.model = model
    return scheme

def burgers(space, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    vspace = space[0]
    pspace = space[1]
    vstorage = storageToSolver(vspace.storage)
    pstorage = storageToSolver(pspace.storage)
    if not (vstorage == pstorage):
        raise KeyError("storages provided differ")

    includes = [ "navierstokes/burgers.cc" ] + vspace._module._includes + pspace._module._includes
    vspaceType = vspace._module._typeName
    pspaceType = pspace._module._typeName
    typeName = "BurgersSchemeWrapper<PRPScheme< " + vspaceType + ", " + pspaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + vspaceType + "::GridPartType, " +\
          vspaceType + "::dimRange+1 " +\
          vspaceType + "::dimRange+1 " +\
        "> > >"

    return module(includes, typeName).Scheme((vspace, pspace), model, name, viscosity, timestep) # ,**kwargs)

from dune.fem.scheme.dgmodel import transform
def dg(model, space=None, penalty=1, solver=None, parameters={},
       penaltyClass=None):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """
    if hasattr(model,"interpolate"):
        warnings.warn("""
        note: the parameter order for the 'schemes' has changes.
              First argument is now the ufl form and the second argument is
              the space.""")
        model,space = space,model

    modelParam = []
    if isinstance(model, (list, tuple)):
        modelParam = model[1:]
        model = model[0]
    if isinstance(model,Equation):
        if space == None:
            try:
                space = model.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no space provided and could not deduce from form provided")
        from dune.fem.model._models import elliptic
        if modelParam:
            model = elliptic(space.grid,model,*modelParam,
                      modelPatch=transform(space,penalty))
        else:
            model = elliptic(space.grid,model,
                      modelPatch=transform(space,penalty))

    spaceType = space._typeName
    useDirichletBC = "true" if model.hasDirichletBoundary else "false"
    modelParam = None
    if isinstance(model, (list, tuple)):
        modelParam = model[1:]
        model = model[0]
    if isinstance(model,Equation):
        from dune.fem.model._models import elliptic
        if modelParam:
            model = elliptic(space.grid,model,*modelParam)
        else:
            model = elliptic(space.grid,model)
    if penaltyClass is None:
        penaltyClass = "DefaultPenalty<"+spaceType+">"
    includes = ["dune/fem/schemes/dgelliptic.hh"]
    operator = lambda linOp,model: "DifferentiableDGEllipticOperator< " +\
                                   ",".join([linOp,model,penaltyClass]) + ">"
    parameters["penalty"] = parameters.get("penalty",penalty)

    includes += ["dune/fem/schemes/diffusionmodel.hh"]
    modelType = "DGDiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters, modelType=modelType)

def dgGalerkin(space, model, penalty, solver=None, parameters={}):
    includes = ["dune/fem/schemes/galerkin.hh"]

    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableDGGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DGDiffusionModelIntegrands<"+model+">"]) + ">"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters)


def _galerkin(integrands, space=None, solver=None, parameters={},
              errorMeasure=None, virtualize=None, schemeName=None ):
    if schemeName is None:
        raise Exception("_galerkin needs a scheme Name: GalerkinScheme or MethodOfLinesScheme")

    if hasattr(integrands,"interpolate"):
        warnings.warn("""
        note: the parameter order for the 'schemes' has changes.
              First argument is now the ufl form and the second argument is
              the space.""")
        integrands,space = space,integrands
    integrandsParam = None
    if isinstance(integrands, (list, tuple)):
        integrandsParam = integrands[1:]
        integrands = integrands[0]
    if isinstance(integrands,Equation):
        if space is None:
            try:
                space = integrands.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no space provided and could not deduce from form provided")
        else:
            try:
                eqSpace = integrands.lhs.arguments()[0].ufl_function_space()
                if not eqSpace.dimRange == space.dimRange:
                    raise ValueError("""range of space used for arguments in equation
                    and the range of space passed to the scheme have to match -
                    space argument to scheme can be 'None'""")
            except AttributeError:
                pass
        from dune.fem.model._models import integrands as makeIntegrands
        if integrandsParam:
            integrands = makeIntegrands(space.grid,integrands,*integrandsParam)
        else:
            integrands = makeIntegrands(space.grid,integrands)
    elif not integrands._typeName.startswith("Integrands"):
        raise ValueError("integrands parameter is not a ufl equation of a integrands model instance")
    if not hasattr(space,"interpolate"):
        raise ValueError("wrong space given")
    from . import module

    storageStr, dfIncludes, dfTypeName, linearOperatorType, defaultSolver,backend = space.storage
    _, solverIncludes, solverTypeName, param = getSolver(solver, space.storage, defaultSolver)

    if virtualize is None:
        virtualize = integrands.virtualized

    includes = [] # integrands._includes
    includes += space._includes + dfIncludes + solverIncludes
    includes += ["dune/fempy/parameter.hh"]
    # molgalerkin includes galerkin.hh so it works for both
    includes += ["dune/fem/schemes/molgalerkin.hh","dune/fem/schemes/dirichletwrapper.hh"]

    spaceType = space._typeName
    valueType = 'std::tuple< typename ' + spaceType + '::RangeType, typename ' + spaceType + '::JacobianRangeType >'
    if virtualize:
        integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + spaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'
    else:
        includes += integrands._includes
        integrandsType = integrands._typeName

    useDirichletBC = "true" if integrands.hasDirichletBoundary else "false"
    typeName = 'Dune::Fem::'+schemeName+'< ' + integrandsType + ', ' +\
            linearOperatorType + ', ' + solverTypeName + ', ' + useDirichletBC + ' >'

    ctors = []
    ctors.append(Constructor(['const ' + spaceType + ' &space', integrandsType + ' &integrands'],
                             ['return new ' + typeName + '( space, std::ref( integrands ) );'],
                             ['"space"_a', '"integrands"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()']))
    ctors.append(Constructor(['const ' + spaceType + ' &space', integrandsType + ' &integrands', 'const pybind11::dict &parameters'],
                             ['return new ' + typeName + '( space, std::ref( integrands ), Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                             ['"space"_a', '"integrands"_a', '"parameters"_a',
                                 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()', 'pybind11::keep_alive< 1, 4 >()']))

    parameters.update(param)
    # scheme = module(includes, typeName, *ctors, backend=backend).Scheme(space, integrands, parameters)
    scheme = module(includes, typeName, backend=backend).Scheme(space, integrands, parameters)
    scheme.model = integrands
    if not errorMeasure is None:
        scheme.setErrorMeasure( errorMeasure );
    return scheme

def galerkin(integrands, space=None, solver=None, parameters={},
             errorMeasure=None, virtualize=None):
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, schemeName='GalerkinScheme')

def molGalerkin(integrands, space=None, solver=None, parameters={},
                errorMeasure=None, virtualize=None):
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, schemeName='MethodOfLinesScheme')


def h1(model, space=None, solver=None, parameters={}):
    """create a scheme for solving second order pdes with continuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """
    if hasattr(model,"interpolate"):
        warnings.warn("""
        note: the parameter order for the 'schemes' has changes.
              First argument is now the ufl form and the second argument is
              the space.""")
        model,space = space,model
    modelParam = None
    if isinstance(model, (list, tuple)):
        modelParam = model[1:]
        model = model[0]
    if isinstance(model,Equation):
        if space == None:
            try:
                space = model.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no space provided and could not deduce from form provided")
        from dune.fem.model._models import elliptic
        if modelParam:
            model = elliptic(space.grid,model,*modelParam)
        else:
            model = elliptic(space.grid,model)

    if not hasattr(space,"interpolate"):
        raise ValueError("wrong space given")

    includes = ["dune/fem/schemes/dirichletwrapper.hh","dune/fem/schemes/elliptic.hh"]

    op = lambda linOp,model: "DifferentiableEllipticOperator< " + ",".join([linOp,model]) + ">"
    if model.hasDirichletBoundary:
        operator = lambda linOp,model: "DirichletWrapperOperator< " + op(linOp,model) + " >"
    else:
        operator = op
    return femschemeModule(space,model,includes,solver,operator,parameters=parameters)

def h1Galerkin(space, model, solver=None, parameters={}):
    from . import module

    includes = [ "dune/fem/schemes/galerkin.hh" ]
    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DiffusionModelIntegrands<"+model+">"]) + ">"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters)


def linearized(scheme, ubar=None, parameters={}):
    from . import module
    schemeType = scheme._typeName
    typeName = "Dune::Fem::LinearizedScheme< " + ", ".join([schemeType]) + " >"
    includes = ["dune/fem/schemes/linearized.hh", "dune/fempy/parameter.hh"] + scheme._includes

    constructor1 = Constructor(['typename DuneType::SchemeType &scheme', 'typename DuneType::DiscreteFunctionType &ubar', 'const pybind11::dict &parameters'],
                               ['return new DuneType( scheme, ubar, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                               ['"scheme"_a', '"ubar"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()'])
    constructor2 = Constructor(['typename DuneType::SchemeType &scheme', 'const pybind11::dict &parameters'],
                               ['return new DuneType( scheme,  Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                               ['"scheme"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()'])
    setup = Method('setup', '&DuneType::setup')

    m = module(includes, typeName, constructor1, constructor2, setup)
    if ubar:
        return m.Scheme(scheme, ubar, parameters)
    else:
        return m.Scheme(scheme, parameters)

def stokes(space, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    vspace = space[0]
    pspace = space[1]
    vstorage = storageToSolver(vspace.storage)
    pstorage = storageToSolver(pspace.storage)
    if not (vstorage == pstorage):
        raise KeyError("storages provided differ")

    includes = [ "navierstokes/stokes.cc" ] + vspace._module._includes + pspace._module._includes
    vspaceType = vspace._module._typeName
    pspaceType = pspace._module._typeName
    typeName = "StokesSchemeWrapper<UzawaScheme< " + vspaceType + ", " + pspaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + vspaceType + "::GridPartType, " +\
          vspaceType + "::dimRange+1 " +\
          vspaceType + "::dimRange+1 " +\
        "> > >"

    return module(includes, typeName).Scheme((vspace, pspace), model, name, viscosity, timestep) #**kwargs)
