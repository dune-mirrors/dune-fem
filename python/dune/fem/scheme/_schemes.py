from __future__ import absolute_import, division, print_function, unicode_literals

import logging
import warnings

from ufl.equation import Equation
from ufl import Form

from dune.generator import Constructor, Method

from dune.fem.operator import _opDirichletIndices

logger = logging.getLogger(__name__)


def _inColor( s, color=None ):
    if color is None:
        return s
    endcolor = '\033[0m'
    return f"{color}{s}{endcolor}"


def _checkNewtonInParameters( params ):
    """
    Args: params a dictionary with parameters

    Returns: a dictionary where every entry containing 'newton' is replaced with
    the same value for a key containing 'nonlinear', unless the same key existed
    before.
    """

    inRed = lambda s :  _inColor(s, '\033[91m')
    inGreen = lambda s : _inColor(s, '\033[92m')

    # check for deprecated use of 'newton.linear' or 'nonlinear.linear'
    warned_linear = False
    parameters = params.copy()
    for key,value in params.items():
        # skip forcing here, warning will appear in NewtonInverseOperator
        if key.find('linear.tolerance.strategy') != -1:
            continue

        if key.find('newton.linear') != -1 or key.find('nonlinear.linear') != -1:
            keypref = ''
            if key.find('newton.') != -1:
                keypref = 'newton.'
                if not warned_linear:
                    warnings.warn(f"""Warning: the parameter key '{inRed("newton.linear")}' is deprecated. Simply remove '{inRed("newton.")}' to avoid this warning!""")
            if key.find('nonlinear.') != -1:
                keypref = 'nonlinear.'
                if not warned_linear:
                    warnings.warn(f"""Warning: the parameter key '{inRed("nonlinear.linear")}' is deprecated. Simply remove '{inRed("nonlinear.")}' to avoid this warning!""")
            # issue warning if this combination has been found
            warned_linear = True
            assert keypref != ''

            newkey = key.replace(keypref, '')
            if newkey in parameters and newkey != "linear.verbose":
                # parameters.pop(key)
                raise KeyError(f"""Parameter dict contains keys '{newkey}' and '{key}'. Mixing new and old parameter keys is not allowed""")
            else:
                parameters[newkey] = parameters.pop(key)

    # check again now for newton only
    warned = False
    params = parameters.copy()
    for key,value in params.items():
        if key.find('newton') != -1:
            # issue warning if newton has been found
            if not warned:
                warnings.warn(f"""Warning: the parameter key '{inRed("newton")}' is deprecated. Replace with '{inGreen("nonlinear")}' to avoid this warning!""")
                warned = True

            newkey = key.replace('newton', 'nonlinear')

            if newkey in parameters:
                parameters.pop(key)
            else:
                parameters[newkey] = parameters.pop(key)

    return parameters

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

    # get storage of solver, it could differ from storage of space
    solverStorage, solver = getSolverStorage(space, solver)

    _, solverIncludes, solverTypeName,param = getSolver(solver,solverStorage, solverStorage.solver)

    includes += ["dune/fem/schemes/femscheme.hh"] +\
                space.cppIncludes + solverStorage.includes + solverIncludes +\
                ["dune/fempy/parameter.hh"]
    spaceType = space.cppTypeName
    if modelType is None:
        includes += ["dune/fem/schemes/diffusionmodel.hh"]
        modelType = "ConservationLawModel< " +\
              "typename " + spaceType + "::GridPartType, " +\
              spaceType + "::dimRange, " +\
              spaceType + "::dimRange, " +\
              "typename " + spaceType + "::RangeFieldType >"
    operatorType = operator(solverStorage.linopType, modelType)
    typeName = "Dune::Fem::FemScheme< " + operatorType + ", " + solverTypeName + " >"
    return includes, typeName

def _linearized(scheme, ubar=None, assemble=True, parameters={}, onlyLinear=True):

    # check for newton in parameters
    parameters = _checkNewtonInParameters( parameters )

    if assemble and not ubar:
        ubar = scheme.space.zero
    if ubar:
        assemble=True
    from . import module
    schemeType = scheme.cppTypeName
    if onlyLinear:
        typeName = "Dune::Fem::LinearScheme< " + ", ".join([schemeType]) + " >"
        baseClasses = ["DuneType::BaseType"]
    else:
        typeName = "Dune::Fem::LinearizedScheme< " + ", ".join([schemeType]) + " >"
        baseClasses = []
    includes = [
        "dune/fem/schemes/linearized.hh",
        "dune/fempy/parameter.hh",
    ] + scheme.cppIncludes

    constructor1 = Constructor(['typename DuneType::SchemeType &scheme', 'typename DuneType::DiscreteFunctionType &ubar', 'const pybind11::dict &parameters'],
                               ['return new DuneType( scheme, ubar, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ) );'],
                               ['"scheme"_a', '"ubar"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()'])
    constructor2 = Constructor(['typename DuneType::SchemeType &scheme', 'const pybind11::dict &parameters'],
                               ['return new DuneType( scheme, Dune::FemPy::pyParameter( "fem.solver.", parameters, std::make_shared< std::string >() ) );'],
                               ['"scheme"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()'])
    if not onlyLinear:
        setup1 = Method('setup', 'static_cast<void (DuneType::*)(const typename DuneType::DiscreteFunctionType &)>(&DuneType::setup)')
        setup2 = Method('setup', 'static_cast<void (DuneType::*)()>(&DuneType::setup)')
    else:
        setup1, setup2 = None,None

    m = module(includes, typeName, constructor1, constructor2, setup1, setup2, baseClasses=baseClasses)

    if onlyLinear:
        from dune.fem.space import addBackend
        # this might differ from space storage
        try:
            backend = scheme._solverBackend
        except AttributeError:
            # if backend has not been added use backend from space
            backend = scheme.space.storage.backend
        if hasattr(m.Scheme,"_backend") and backend is not None:
            if backend == 'as_numpy':
                from scipy.sparse import csr_matrix
            addBackend(m.Scheme,backend)

    # Remove a Prefix in Python
    rmPre = lambda txt,pre: txt if not txt.startswith(pre) else txt[len(pre):]
    p = dict( [[rmPre(k,"linear."),v] for k,v in scheme.parameters.items()] )
    params = { **p, **parameters }

    if assemble:
        linearizedScheme = m.Scheme(scheme, ubar=ubar, parameters=params)
    else:
        linearizedScheme = m.Scheme(scheme, parameters=params)

    linearizedScheme.nonLinearModel = scheme.model
    linearizedScheme.parameters = params
    linearizedScheme.__class__.dirichletIndices = _opDirichletIndices

    return linearizedScheme

def linearized(scheme, ubar=None, assemble=True, parameters={}):
    return _linearized(scheme,ubar,assemble,parameters,onlyLinear=False)

def femschemeModule(space, model, includes, solver, operator, *args,
        parameters={},
        modelType = None, ctorArgs={}):

    # check for newton in parameters
    parameters = _checkNewtonInParameters( parameters )

    from . import module
    # get storage of solver, it could differ from storage of space
    solverStorage, solver = getSolverStorage(space, solver)
    *_, param = getSolver(solver, solverStorage, solverStorage.solver)

    includes, typeName = femscheme(includes, space, solver, operator, modelType)
    parameters.update(param)
    mod = module(includes, typeName, *args, backend=solverStorage.backend)
    scheme = mod.Scheme(space, model, parameters=parameters, **ctorArgs)
    scheme.model = model
    scheme.parameters = parameters
    scheme._solverBackend = solverStorage.backend
    scheme.__class__.linear = _linearized
    scheme.__class__.dirichletIndices = _opDirichletIndices
    return scheme

def dg(model, space=None, penalty=1, solver=None, parameters={},
       penaltyClass=None):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from dune.fem.scheme.dgmodel import transform
    # check for newton in parameters
    parameters = _checkNewtonInParameters( parameters )

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

def _massLumpingGalerkin(integrands, integrandsParam=None, massIntegrands=None, space=None, solver=None, parameters={},
                         errorMeasure=None, virtualize=None, **kwargs):
    """
    Parameters:

        integrands     Model of the weak form of the PDEs solved by Galerkin scheme
        massIntegrands Model of the mass terms to be lumped
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
    assert massIntegrands is not None, "Missing mass term for massLumpingGalerkin"
    _schemeName = "MassLumpingScheme"

    # check for newton in parameters
    parameters = _checkNewtonInParameters( parameters )

    if isinstance(integrands,Equation):
        if space is None:
            try:
                space = integrands.lhs.arguments()[0].ufl_function_space()
            except AttributeError:
                raise ValueError("no space provided and could not deduce from form provided")
            except IndexError:
                from dune.ufl import Constant
                from ufl import TrialFunction, TestFunction, dx, inner
                # the following is a hack to avoid issues with a missing 'normal' integrand
                space = massIntegrands.lhs.arguments()[0].ufl_function_space()
                integrands = Constant(0)*inner( TrialFunction(space), TestFunction(space) ) * dx == 0
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
    else:
        try:
            if not integrands.cppTypeName.startswith("Integrands"):
                raise ValueError("integrands parameter is not a ufl equation of a integrands model instance")
        except AttributeError:
            raise ValueError("first argument should be a ufl equation (not only a form) or an 'integrands' model")
    if not hasattr(space,"interpolate"):
        raise ValueError("wrong space given")

    if isinstance(massIntegrands, Equation):
        from dune.fem.model._models import integrands as makeIntegrands
        massIntegrands = makeIntegrands(space.gridView, massIntegrands)

    from . import module

    # get storage of solver, it could differ from storage of space
    solverStorage, solver = getSolverStorage(space, solver)

    _, solverIncludes, solverTypeName, param = getSolver(solver, solverStorage, solverStorage.solver)

    # check if parameters have an option preconditioning and if this is a callable
    preconditioning = None
    precondkey = 'linear.preconditioning.method'
    if precondkey in parameters:
        # if preconditioning is callable then store as preconditioning
        # and remove from parameters
        if callable(parameters[ precondkey ]):
            # store as preconditioning object and remove from dict
            preconditioning = parameters[ precondkey ]
            parameters.pop( precondkey )

    assert integrands.virtualized and massIntegrands.virtualized, "MassLumping only works with virtualized integrands!"
    virtualize = True

    includes = [] # integrands.cppIncludes
    includes += space.cppIncludes + space.storage.includes + solverIncludes
    includes += ["dune/fempy/parameter.hh"]
    # molgalerkin includes galerkin.hh so it works for both
    includes += ["dune/fem/schemes/masslumping.hh","dune/fem/schemes/dirichletwrapper.hh"]

    spaceType = space.cppTypeName
    valueType = 'std::tuple< typename ' + spaceType + '::RangeType, typename ' + spaceType + '::JacobianRangeType >'
    if virtualize:
        integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + spaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'
        massIntegrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + spaceType + '::GridPartType, ' + massIntegrands._domainValueType + ", " + massIntegrands._rangeValueType+ ' >'
    else:
        includes += integrands.cppIncludes
        includes += massIntegrands.cppIncludes
        integrandsType = integrands.cppTypeName
        massIntegrandsType = massIntegrands.cppTypeName

    useDirichletBC = "true" if integrands.hasDirichletBoundary else "false"
    typeName = 'Dune::Fem::'+_schemeName+'< ' + integrandsType + ', ' + massIntegrandsType + ', ' +\
            solverStorage.linopType + ', ' + solverTypeName + ', ' + useDirichletBC + ' >'

    parameters.update(param)
    scheme = module(includes, typeName, backend=solverStorage.backend).Scheme(space, integrands, massIntegrands, parameters)
    scheme.model = integrands
    scheme.massModel = massIntegrands
    # store solver backend
    scheme._solverBackend = solverStorage.backend
    scheme.parameters = parameters

    scheme.__class__.linear = _linearized
    scheme.__class__.dirichletIndices = _opDirichletIndices

    if not errorMeasure is None:
        scheme.setErrorMeasure( errorMeasure );

    # if preconditioning was passed as callable then store in scheme, otherwise None is stored
    scheme.preconditioning = preconditioning

    try:
        from dune.fem import parameter
        logTag = parameters["logging"]
        scheme.parameterLog = lambda: parameter.log()[logTag]
    except KeyError:
        pass

    return scheme

def _galerkin(integrands, space=None, solver=None, parameters={},
              errorMeasure=None, virtualize=None, _schemeName=None, **kwargs ):
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

    # check for newton in parameters
    parameters = _checkNewtonInParameters( parameters )

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

    ################################################################
    ##
    ## Scan metadata of integrands to detect which scheme to use
    ##
    ################################################################
    mass = None
    other = []
    if isinstance(integrands,Equation):
        # move entire equation to lhs to have a consistent representation
        form = integrands.lhs - integrands.rhs
        s = set( [ tuple(i.metadata().items()) for i in form.integrals() ] )
        formVec = {}
        lumped = ("quadrature_rule", "lumped")
        for i in form.integrals():
            formVec[tuple(i.metadata().items())] = formVec.get(tuple(i.metadata().items()),[]) + [i]

        for dx in s:
            if len(dx) > 0:
                if dx[0] == (lumped):
                    mass = Form(formVec[dx])
                else:
                    raise ValueError("currently only quadrature_rule implemented is 'lumped'")
            else:
                other += formVec[dx]
        other = Form(other)

    # if mass was given use _massLumpingGalerkin else continue
    if mass is not None:
        mass = mass == 0
        integrands = other == 0
        return _massLumpingGalerkin(integrands, integrandsParam=integrandsParam,
                                    massIntegrands=mass, space=space, solver=solver,
                                    parameters=parameters,
                                    errorMeasure=errorMeasure,
                                    virtualize=virtualize)

    ## GalerkinOperator with one integrand

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
    else:
        try:
            if not integrands.cppTypeName.startswith("Integrands"):
                raise ValueError("integrands parameter is not a ufl equation of a integrands model instance")
        except AttributeError:
            raise ValueError("first argument should be a ufl equation (not only a form) or an 'integrands' model")
    if not hasattr(space,"interpolate"):
        raise ValueError("wrong space given")
    from . import module

    # get storage of solver, it could differ from storage of space
    solverStorage, solver = getSolverStorage(space, solver)

    _, solverIncludes, solverTypeName, param = getSolver(solver, solverStorage, solverStorage.solver)

    # check if parameters have an option preconditioning and if this is a callable
    preconditioning = None
    precondkey = 'linear.preconditioning.method'
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
    includes += space.cppIncludes + space.storage.includes + solverIncludes
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
            solverStorage.linopType + ', ' + solverTypeName + ', ' + useDirichletBC + ' >'

    parameters.update(param)
    scheme = module(includes, typeName, backend=solverStorage.backend).Scheme(space, integrands, parameters)
    scheme.model = integrands
    if not errorMeasure is None:
        scheme.setErrorMeasure( errorMeasure );

    # if preconditioning was passed as callable then store in scheme, otherwise None is stored
    scheme.preconditioning = preconditioning
    # store solver backend
    scheme._solverBackend = solverStorage.backend

    scheme.parameters = parameters

    scheme.__class__.linear = _linearized
    scheme.__class__.dirichletIndices = _opDirichletIndices

    try:
        from dune.fem import parameter
        logTag = parameters["logging"]
        scheme.parameterLog = lambda: parameter.log()[logTag]
    except KeyError:
        pass

    return scheme

def galerkin(integrands, space=None, solver=None, parameters={},
             errorMeasure=None, virtualize=None, **kwargs):
    domainSpace = kwargs.get("domainSpace",space)
    rangeSpace = kwargs.get("rangeSpace",space)
    if space is None:
        space = domainSpace
    if domainSpace != rangeSpace or domainSpace != space:
        raise ValueError("provided spaces not compatible")

    galerkin.__doc__ = _galerkin.__doc__
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, _schemeName='GalerkinScheme', **kwargs)

def molGalerkin(integrands, space=None, solver=None, parameters={},
                errorMeasure=None, virtualize=None, **kwargs):
    molGalerkin.__doc__ = _galerkin.__doc__
    return _galerkin(integrands, space=space, solver=solver,
                     parameters=parameters, errorMeasure=errorMeasure,
                     virtualize=virtualize, _schemeName='MethodOfLinesScheme', **kwargs)


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
