from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

def burgers(space_or_df, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, spaceAndStorage
    storage = kwargs.pop("storage","fem")
    vspace, vstorage = spaceAndStorage(space_or_df[0],storage)
    pspace, pstorage = spaceAndStorage(space_or_df[1],storage)
    if not (vstorage == pstorage):
        raise KeyError("storages provided differe")
    storage = vstorage

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError(\
            "Parameter error in FemScheme with "+\
            "storage=" + storage)

    includes = [ "navierstokes/burgers.cc" ] + vspace._module._includes + pspace._module._includes
    vspaceType = vspace._module._typeName
    pspaceType = pspace._module._typeName
    typeName = "BurgersSchemeWrapper<PRPScheme< " + vspaceType + ", " + pspaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + vspaceType + "::GridPartType, " +\
          vspaceType + "::dimRange+1 " +\
        "> > >"

    return module(storage, includes, typeName).Scheme((vspace,pspace),model,name,viscosity,timestep) # ,**kwargs)

def dg(space_or_df, model, name="tmp", **kwargs):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, spaceAndStorage
    storage = kwargs.pop("storage","fem")
    space, storage = spaceAndStorage(space_or_df,storage)

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError(\
            "Parameter error in FemScheme with "+\
            "storage=" + storage)

    includes = [ "dune/fem/schemes/dgelliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableDGEllipticOperator, " +\
          storage + " >"

    return module(storage, includes, typeName).Scheme(space,model,name,**kwargs)


def galerkin(space, model, name="tmp", storage="fem", parameters={}):
    from . import module, spaceAndStorage
    space, storage = spaceAndStorage(space, storage)

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Eigen" or storage == "eigen":
        storage = "eigen"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError("Parameter error in Galerkin scheme with storage=" + storage)

    spaceType = space._module._typeName
    gridPartType = "typename " + spaceType + "::GridPartType"
    dimRange = spaceType + "::dimRange"
    rangeFieldType = "typename " + spaceType + "::RangeFieldType"
    modelType = "IntegrandsModel< " + ", ".join([gridPartType, dimRange, rangeFieldType]) + " >"

    typeName = "Dune::Fem::GalerkinScheme< " + ", ".join([spaceType, modelType, storage]) + " >"
    includes = ["dune/fem/schemes/galerkin.hh"] + space._module._includes

    return module(storage, includes, typeName).Scheme(space, model, name, parameters)


def h1(space, model, storage="fem", parameters={}):
    """create a scheme for solving second order pdes with continuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, spaceAndStorage
    space, storage = spaceAndStorage(space,storage)

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Eigen" or storage == "eigen":
        try:
            checkconfiguration.preprocessorTest([ ("#if HAVE_EIGEN","Eigen package is not available") ])
        except checkconfiguration.ConfigurationError as err:
            if logger.getEffectiveLevel() == logging.DEBUG:
                raise
            else:
                print("configuration error while creating a discrete function with storage=eigen exiting...")
                sys.exit(-1)
        storage = "eigen"

    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError(\
            "Parameter error in FemScheme with "+\
            "storage=" + storage)

    includes = [ "dune/fem/schemes/elliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableEllipticOperator, " +\
          storage + " >"

    return module(storage, includes, typeName).Scheme(space,model,parameters)


def h1Galerkin(space, model, storage="fem", parameters={}):
    from . import module, spaceAndStorage
    space, storage = spaceAndStorage(space, storage)

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Eigen" or storage == "eigen":
        storage = "eigen"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError("Parameter error in Galerkin scheme with storage=" + storage)

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

    return module(storage, includes, typeName, [constructor]).Scheme(space, model, parameters)


def nvdg(space_or_df, model, name="tmp", **kwargs):
    """create a scheme for solving non variational second order pdes with discontinuous finite element

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, spaceAndStorage
    storage = kwargs.pop("storage","fem")
    space, storage = spaceAndStorage(space_or_df,storage)

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError(\
            "Parameter error in FemScheme with "+\
            "storage=" + storage)

    includes = [ "dune/fem/schemes/nvdgelliptic.hh", "dune/fem/schemes/femscheme.hh" ] + space._module._includes
    spaceType = space._module._typeName
    typeName = "FemScheme< " + spaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + spaceType + "::GridPartType, " +\
          spaceType + "::dimRange, " +\
          "typename " + spaceType + "::RangeFieldType >, DifferentiableNVDGEllipticOperator, " +\
          storage + " >"

    return module(storage, includes, typeName).Scheme(space,model,name,**kwargs)

def stokes(space_or_df, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """

    from . import module, spaceAndStorage
    storage = kwargs.pop("storage","fem")
    vspace, vstorage = spaceAndStorage(space_or_df[0],storage)
    pspace, pstorage = spaceAndStorage(space_or_df[1],storage)
    if not (vstorage == pstorage):
        raise KeyError("storages provided differe")
    storage = vstorage

    if storage == "Adaptive" or storage == "adaptive":
        storage = "fem"
    elif storage == "Istl" or storage == "istl":
        storage = "istl"
    elif storage == "Numpy" or storage == "numpy":
        storage = "numpy"
    elif storage == "Fem" or storage == "fem":
        storage = "fem"
    else:
        raise KeyError(\
            "Parameter error in FemScheme with "+\
            "storage=" + storage)

    includes = [ "navierstokes/stokes.cc" ] + vspace._module._includes + pspace._module._includes
    vspaceType = vspace._module._typeName
    pspaceType = pspace._module._typeName
    typeName = "StokesSchemeWrapper<UzawaScheme< " + vspaceType + ", " + pspaceType + ", " +\
        "DiffusionModel< " +\
          "typename " + vspaceType + "::GridPartType, " +\
          vspaceType + "::dimRange+1 " +\
        "> > >"

    return module(storage, includes, typeName).Scheme((vspace,pspace),model,name,viscosity,timestep) #**kwargs)
