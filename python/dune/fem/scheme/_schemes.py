from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
import warnings

from ufl.equation import Equation

from dune.generator import Constructor, Method
from dune.common.utility import isString

logger = logging.getLogger(__name__)

def getSolverStorage(space, solver):
    """
    Return storage that fits chosen solver (default: space.storage)
    """
    from dune.fem.discretefunction import _storage
    storage = space.storage

    # this feature only works for numpy storage of the space
    if storage[0] == 'numpy' and isinstance(solver, tuple) and len(solver) > 1:
        try:
            # changing the storage only works for petsc and istl, but is not needed for suitesparse
            return _storage(dfStorage=storage[0], solverStorage=solver[0])(space), solver[1]
        except ValueError:
            pass

    # otherwise simply return arguments that were passed
    return storage, solver

def getSolver(solver, storage, default):
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
                space.cppIncludes + dfIncludes + solverIncludes +\
                ["dune/fempy/parameter.hh"]
    spaceType = space.cppTypeName
    if modelType is None:
        includes += ["dune/fem/schemes/diffusionmodel.hh"]
        modelType = "ConservationLawModel< " +\
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
            model = elliptic(space.gridView,model,*modelParam,
                      modelPatch=transform(space,penalty))
        else:
            model = elliptic(space.gridView,model,
                      modelPatch=transform(space,penalty))

    spaceType = space.cppTypeName
    useDirichletBC = "true" if model.hasDirichletBoundary else "false"
    modelParam = None
    if isinstance(model, (list, tuple)):
        modelParam = model[1:]
        model = model[0]
    if isinstance(model,Equation):
        from dune.fem.model._models import elliptic
        if modelParam:
            model = elliptic(space.gridView,model,*modelParam)
        else:
            model = elliptic(space.gridView,model)
    if penaltyClass is None:
        penaltyClass = "DefaultPenalty<"+spaceType+">"
    includes = ["dune/fem/schemes/dgelliptic.hh"]
    operator = lambda linOp,model: "DifferentiableDGEllipticOperator< " +\
                                   ",".join([linOp,model,penaltyClass]) + ">"
    parameters["penalty"] = parameters.get("penalty",penalty)

    includes += ["dune/fem/schemes/conservationlawmodel.hh"]
    modelType = "DGConservationLawModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters, modelType=modelType)

def dgGalerkin(space, model, penalty, solver=None, parameters={}):
    includes = ["dune/fem/schemes/galerkin.hh"]

    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableDGGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DGConservationLawModelIntegrands<"+model+">"]) + ">"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters)


def _galerkin(integrands, space=None, solver=None, parameters={},
              errorMeasure=None, virtualize=None, _schemeName=None ):
    """
    Parameters:

        integrands   Model of the weak form of the PDEs solved by Galerkin scheme
        space        Discrete function space
        solver       String (e.g. 'gmres', 'bicgstab', 'cg'...) or tuple
                     containing backend and solver name,
                     e.g. ('petsc', 'gmres') or ('istl', 'bicgstab') or
                     ('suitesparse', 'umfpack'). See documentation for complete
                     list.
        parameters   dictionary with parameters passed to the nonlinear and linear solvers
        errorMeasure Overloading the nonlinear solver stopping criterion,
                     i.e. a function `f( w, dw, float )` where w and dw are discrete functions returning a bool whether the
                     solver has converged or not.
        virtualize   If True, integrands will be virtualized to avoid
                     re-compilation. (default: True)
    """

    if _schemeName is None:
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
            integrands = makeIntegrands(space.gridView,integrands,*integrandsParam)
        else:
            integrands = makeIntegrands(space.gridView,integrands)
    elif not integrands.cppTypeName.startswith("Integrands"):
        raise ValueError("integrands parameter is not a ufl equation of a integrands model instance")
    if not hasattr(space,"interpolate"):
        raise ValueError("wrong space given")
    from . import module

    storageStr, dfIncludes, dfTypeName, _, _, _ = space.storage

    # get storage of solver, it could differ from storage of space
    solverStorage, solver = getSolverStorage(space, solver)

    _, _, _, linearOperatorType, defaultSolver, backend = solverStorage
    _, solverIncludes, solverTypeName, param = getSolver(solver, solverStorage, defaultSolver)

    # check if parameters have an option preconditioning and if this is a callable
    preconditioning = None
    precondkey = 'newton.linear.preconditioning.method'
    if precondkey in parameters:
        # if preconditioning is callable then store as preconditioning
        # and remove from parameters
        if callable(parameters[ precondkey ]):
            # store as preconditioning object and remove from dict
            preconditioning = parameters[ precondkey ]
            parameters.pop( precondkey )

    if virtualize is None:
        virtualize = integrands.virtualized

    includes = [] # integrands.cppIncludes
    includes += space.cppIncludes + dfIncludes + solverIncludes
    includes += ["dune/fempy/parameter.hh"]
    # molgalerkin includes galerkin.hh so it works for both
    includes += ["dune/fem/schemes/molgalerkin.hh","dune/fem/schemes/dirichletwrapper.hh"]

    spaceType = space.cppTypeName
    valueType = 'std::tuple< typename ' + spaceType + '::RangeType, typename ' + spaceType + '::JacobianRangeType >'
    if virtualize:
        integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + spaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'
    else:
        includes += integrands.cppIncludes
        integrandsType = integrands.cppTypeName

    useDirichletBC = "true" if integrands.hasDirichletBoundary else "false"
    typeName = 'Dune::Fem::'+_schemeName+'< ' + integrandsType + ', ' +\
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

    # if preconditioning was passed as callable then store in scheme, otherwise None is stored
    scheme.preconditioning = preconditioning

    return scheme

def galerkin(integrands, space=None, solver=None, parameters={},
             errorMeasure=None, virtualize=None):
    galerkin.__doc__ = _galerkin.__doc__
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, _schemeName='GalerkinScheme')

def molGalerkin(integrands, space=None, solver=None, parameters={},
                errorMeasure=None, virtualize=None):
    molGalerkin.__doc__ = _galerkin.__doc__
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, _schemeName='MethodOfLinesScheme')


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
            model = elliptic(space.gridView,model,*modelParam)
        else:
            model = elliptic(space.gridView,model)

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
            ",".join([linOp,"Dune::Fem::ConservationLawModelIntegrands<"+model+">"]) + ">"

    return femschemeModule(space,model,includes,solver,operator,parameters=parameters)


def linearized(scheme, ubar=None, parameters={}):
    from . import module
    schemeType = scheme.cppTypeName
    typeName = "Dune::Fem::LinearizedScheme< " + ", ".join([schemeType]) + " >"
    includes = ["dune/fem/schemes/linearized.hh", "dune/fempy/parameter.hh"] + scheme.cppIncludes

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
