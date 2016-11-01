from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

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

def dg(space, model, name="tmp", **kwargs):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    includes = [ "dune/fem/schemes/dgelliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableDGEllipticOperator, " +\
          storage + " >"

    return module(includes, typeName).Scheme(space,model,name,**kwargs)


def dgGalerkin(space, model, penalty, parameters={}):
    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    spaceType = space._module._typeName
    gridPartType = "typename " + spaceType + "::GridPartType"
    dimRange = spaceType + "::dimRange"
    rangeFieldType = "typename " + spaceType + "::RangeFieldType"
    modelType = "DiffusionModel< " + ", ".join([gridPartType, dimRange, rangeFieldType]) + " >"
    integrandsType = "Dune::Fem::DGDiffusionModelIntegrands< " + modelType + " >"

    typeName = "Dune::Fem::GalerkinScheme< " + ", ".join([spaceType, integrandsType, storage]) + " >"
    includes = ["dune/fem/schemes/galerkin.hh", "dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"] + space._module._includes

    constructor = ['[] ( ' + typeName + ' &self, const ' + spaceType + ' &space, const ' + modelType + ' &model, const typename ' + modelType + '::RangeFieldType &penalty, const pybind11::dict &parameters ) {',
                   '    new (&self) ' + typeName + '( space, ' + integrandsType + '( model, penalty ), Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );',
                   '  }, "space"_a, "model"_a, "penalty"_a, "parameters"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >()']

    return module(includes, typeName, [constructor]).Scheme(space, model, penalty, parameters)


def galerkin(space, model, parameters={}):
    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    spaceType = space._module._typeName
    gridPartType = "typename " + spaceType + "::GridPartType"
    dimRange = spaceType + "::dimRange"
    rangeFieldType = "typename " + spaceType + "::RangeFieldType"
    modelType = "IntegrandsModel< " + ", ".join([gridPartType, dimRange, rangeFieldType]) + " >"

    typeName = "Dune::Fem::GalerkinScheme< " + ", ".join([spaceType, modelType, storage]) + " >"
    includes = ["dune/fem/schemes/galerkin.hh"] + space._module._includes

    return module(includes, typeName).Scheme(space, model, parameters)


def h1(space, model, parameters={}):
    """create a scheme for solving second order pdes with continuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    includes = [ "dune/fem/schemes/elliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._includes
    spaceType = space._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableEllipticOperator, " +\
          storage + " >"

    return module(includes, typeName).Scheme(space,model,parameters)


def h1Galerkin(space, model, parameters={}):
    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    spaceType = space._module._typeName
    gridPartType = "typename " + spaceType + "::GridPartType"
    dimRange = spaceType + "::dimRange"
    rangeFieldType = "typename " + spaceType + "::RangeFieldType"
    modelType = "DiffusionModel< " + ", ".join([gridPartType, dimRange, rangeFieldType]) + " >"
    integrandsType = "Dune::Fem::DiffusionModelIntegrands< " + modelType + " >"

    typeName = "Dune::Fem::GalerkinScheme< " + ", ".join([spaceType, integrandsType, storage]) + " >"
    includes = ["dune/fem/schemes/galerkin.hh", "dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"] + space._module._includes

    constructor = ['[] ( ' + typeName + ' &self, const ' + spaceType + ' &space, const ' + modelType + ' &model, const pybind11::dict &parameters ) {',
                   '    new (&self) ' + typeName + '( space, ' + integrandsType + '( model ), Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );',
                   '  }, "space"_a, "model"_a, "parameters"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >()']

    return module(includes, typeName, [constructor]).Scheme(space, model, parameters)


def linearized(scheme, ubar=None, parameters={}):
    from . import module, storageToSolver
    schemeType = scheme._typeName
    typeName = "Dune::Fem::LinearizedScheme< " + ", ".join([schemeType]) + " >"
    includes = ["dune/fem/schemes/linearized.hh", "dune/fempy/parameter.hh"] + scheme._includes

    constructor1 = ['[] ( DuneType &self,',
                         'typename DuneType::SchemeType &scheme,',
                         'typename DuneType::DiscreteFunctionType &ubar,',
                         'const pybind11::dict &parameters ) {',
                   '   new (&self) DuneType( scheme, ubar, Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );',
                   '  }, "scheme"_a, "ubar"_a, "parameters"_a,',
                   '     pybind11::keep_alive< 1, 2 >()']
    constructor2 = ['[] ( DuneType &self,',
                         'typename DuneType::SchemeType &scheme,',
                         'const pybind11::dict &parameters ) {',
                   '   new (&self) DuneType( scheme,  Dune::FemPy::pyParameter( parameters, std::make_shared< std::string >() ) );',
                   '  }, "scheme"_a, "parameters"_a,',
                   '     pybind11::keep_alive< 1, 2 >()']
    method = ['setup', 'DuneType::setup']

    m = module(includes, typeName, [constructor1,constructor2], [method])
    if ubar:
        return m.Scheme(scheme, ubar, parameters)
    else:
        return m.Scheme(scheme, parameters)

def nvdg(space, model, name="tmp", **kwargs):
    """create a scheme for solving non variational second order pdes with discontinuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    storage = storageToSolver(space.storage)

    includes = [ "dune/fem/schemes/nvdgelliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableNVDGEllipticOperator, " +\
          storage + " >"

    return module(includes, typeName).Scheme(space, model, name, **kwargs)


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
