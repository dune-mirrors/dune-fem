from __future__ import absolute_import, division, print_function, unicode_literals

from . import module, spaceAndStorage

def create(space_or_df, model, name, viscosity, timestep, **kwargs):
    """create a scheme for solving quasi stokes type saddle point problem with continuous finite-elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """
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

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
