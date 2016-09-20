from __future__ import absolute_import, division, print_function, unicode_literals

from . import module, spaceAndStorage

def create(space_or_df, model, name="tmp", **kwargs):
    """create a scheme for solving second order pdes with discontinuous finite elements

    Args:

    Returns:
        Scheme: the constructed scheme
    """
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

if __name__ == "__main__":
    import doctest
    doctest.testmod(optionflags=doctest.ELLIPSIS)
