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

def h1(space_or_df, model, name="tmp", **kwargs):
    """create a scheme for solving second order pdes with continuous finite element

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
    elif storage == "Eigen" or storage == "eigen":
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

    return module(storage, includes, typeName).Scheme(space,model,name,**kwargs)

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
