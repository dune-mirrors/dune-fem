from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

def femscheme(includes, space, operator):
    from . import storageToSolver
    import dune.create as create
    storage = storageToSolver(space.storage)
    dfIncludes, dfTypeName, linearOperatorType = create.discretefunction(space.storage,space,instanciate=False)
    includes += ["dune/fem/schemes/femscheme.hh"] +\
                space._module._includes + dfIncludes +\
                ["dune/fem/schemes/diffusionmodel.hh", "dune/fempy/parameter.hh"]
    spaceType = space._module._typeName
    modelType = "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >"
    operatorType = operator(linearOperatorType,modelType)
    typeName = "FemScheme< " + operatorType + ", " + storage + " >"
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
    vstorage = storageToSolver(vspape.storage)
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

    return module(includes, typeName).Scheme((vspace,pspace),model,name,viscosity,timestep) # ,**kwargs)

def dg(space, model, **kwargs):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module

    includes = [ "dune/fem/schemes/dgelliptic.hh"]
    operator = lambda linOp,model: "DifferentiableDGEllipticOperator< " +\
                                   ",".join([linOp,model]) + ">"
    includes, typeName = femscheme(includes, space, operator)

    return module(includes, typeName).Scheme(space,model,**kwargs)

def dgGalerkin(space, model, penalty, parameters={}):
    from . import module

    includes = [ "dune/fem/schemes/galerkin.hh" ]

    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableDGGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DGDiffusionModelIntegrands<"+model+">"]) + ">"

    includes, typeName = femscheme(includes, space, operator)

    return module(includes, typeName).Scheme(space, model, parameters)

def h1(space, model, parameters={}):
    """create a scheme for solving second order pdes with continuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """
    from . import module
    includes = [ "dune/fem/schemes/elliptic.hh" ]

    operator = lambda linOp,model: "DifferentiableEllipticOperator< " +\
                                   ",".join([linOp,model]) + ">"
    includes, typeName = femscheme(includes, space, operator)

    return module(includes, typeName).Scheme(space,model,parameters)

def h1Galerkin(space, model, parameters={}):
    from . import module

    includes = [ "dune/fem/schemes/galerkin.hh" ]
    operator = lambda linOp,model: "Dune::Fem::ModelDifferentiableGalerkinOperator< " +\
            ",".join([linOp,"Dune::Fem::DiffusionModelIntegrands<"+model+">"]) + ">"

    includes, typeName = femscheme(includes, space, operator)

    return module(includes, typeName).Scheme(space, model, parameters)


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

    return module(includes, typeName).Scheme(space,model,name,**kwargs)


def stokes(space_or_df, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, storageToSolver
    vspace = space[0]
    pspace = space[1]
    vstorage = storageToSolver(vspape.storage)
    pstorage = storageToSolver(pspace.storage)
    if not (vstorage == pstorage):
        raise KeyError("storages provided differe")

    includes = [ "navierstokes/stokes.cc" ] + vspace._module._includes + pspace._module._includes
    vspaceType = vspace._module._typeName
    pspaceType = pspace._module._typeName
    typeName = "StokesSchemeWrapper<UzawaScheme< " + vspaceType + ", " + pspaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + vspaceType + "::GridPartType, " +\
          vspaceType + "::dimRange+1 " +\
        "> > >"

    return module(includes, typeName).Scheme((vspace,pspace),model,name,viscosity,timestep) #**kwargs)
