from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging

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

def femscheme(includes, space, solver, operator):
    storageStr, dfIncludes, dfTypeName, linearOperatorType, defaultSolver = space.storage
    _, solverIncludes, solverTypeName = getSolver(solver,space.storage,defaultSolver)

    includes += ["dune/fem/schemes/femscheme.hh"] +\
                space._includes + dfIncludes + solverIncludes +\
                ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]
    spaceType = space._typeName
    modelType = "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"
    operatorType = operator(linearOperatorType, modelType)
    typeName = "FemScheme< " + operatorType + ", " + solverTypeName + " >"
    return includes, typeName

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
        "> > >"

    return module(includes, typeName).Scheme((vspace, pspace), model, name, viscosity, timestep) # ,**kwargs)

def dg(space, model, solver=None, **kwargs):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module

    includes = ["dune/fem/schemes/dgelliptic.hh"]
    operator = lambda linOp,model: "DifferentiableDGEllipticOperator< " +\
                                   ",".join([linOp,model]) + ">"
    includes, typeName = femscheme(includes, space, solver, operator)

    return module(includes, typeName).Scheme(space,model,**kwargs)

def dgGalerkin(space, model, penalty, solver=None, parameters={}):
    from . import module

    includes = ["dune/fem/schemes/galerkin.hh"]

    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableDGGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DGDiffusionModelIntegrands<"+model+">"]) + ">"

    includes, typeName = femscheme(includes, space, solver, operator)

    return module(includes, typeName).Scheme(space, model, parameters)


def galerkin(space, integrands, solver=None, parameters={}):
    integrandsParam = None
    if isinstance(integrands, (list, tuple)):
        integrandsParam = integrands[1:]
        integrands = integrands[0]
    if isinstance(integrands,Equation):
        from dune.fem.model._models import integrands as makeIntegrands
        if integrandsParam:
            integrands = makeIntegrands(space.grid,integrands,*integrandsParam)
        else:
            integrands = makeIntegrands(space.grid,integrands)
    from . import module

    storageStr, dfIncludes, dfTypeName, linearOperatorType, defaultSolver = space.storage
    _, solverIncludes, solverTypeName = getSolver(solver, space.storage, defaultSolver)

    includes = ["dune/fem/schemes/galerkin.hh"]
    includes += space._includes + dfIncludes + solverIncludes
    includes += ["dune/fempy/parameter.hh"]

    spaceType = space._typeName
    valueType = 'std::tuple< typename ' + spaceType + '::RangeType, typename ' + spaceType + '::JacobianRangeType >'
    integrandsType = 'Dune::Fem::VirtualizedIntegrands< typename ' + spaceType + '::GridPartType, ' + integrands._domainValueType + ", " + integrands._rangeValueType+ ' >'

    typeName = 'Dune::Fem::GalerkinScheme< ' + integrandsType + ', ' + linearOperatorType + ', ' + solverTypeName + ' >'

    ctors = []
    ctors.append(Constructor(['const ' + spaceType + ' &space', integrandsType + ' &integrands'],
                             ['return new ' + typeName + '( space, std::ref( integrands ) );'],
                             ['"space"_a', '"integrands"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()']))
    ctors.append(Constructor(['const ' + spaceType + ' &space', integrandsType + ' &integrands', 'const pybind11::dict &parameters'],
                             ['return new ' + typeName + '( space, std::ref( integrands ), Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );'],
                             ['"space"_a', '"integrands"_a', '"parameters"_a', 'pybind11::keep_alive< 1, 2 >()', 'pybind11::keep_alive< 1, 3 >()']))

    return module(includes, typeName, *ctors).Scheme(space, integrands, parameters)


def h1(space, model, solver=None, parameters={}):
    """create a scheme for solving second order pdes with continuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """
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

    from . import module
    includes = ["dune/fem/schemes/elliptic.hh"]

    operator = lambda linOp,model: "DifferentiableEllipticOperator< " +\
                                   ",".join([linOp,model]) + ">"
    includes, typeName = femscheme(includes, space, solver, operator)

    return module(includes, typeName).Scheme(space,model,parameters)


def h1Galerkin(space, model, solver=None, parameters={}):
    from . import module

    includes = [ "dune/fem/schemes/galerkin.hh" ]
    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DiffusionModelIntegrands<"+model+">"]) + ">"

    includes, typeName = femscheme(includes, space, solver, operator)

    return module(includes, typeName).Scheme(space, model, parameters)


def linearized(scheme, ubar=None, solver=None, parameters={}):
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
        "> > >"

    return module(includes, typeName).Scheme((vspace, pspace), model, name, viscosity, timestep) #**kwargs)
