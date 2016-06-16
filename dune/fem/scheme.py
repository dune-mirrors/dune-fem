"""Functions for creating python modules and C++ classes for schemes.
"""

from __future__ import print_function
import importlib
import subprocess
import hashlib
import os.path
import re

from .. import femmpi
from ..generator import generator

myGenerator = generator.Generator("Scheme")

def getModule(scheme, **parameters):
    """Create a scheme module using the scheme-database.

    Takes its arguments from get(). This function uses the given arguments to create a python module called
    duneschemexxx.py where xxx is a number derived from the scheme type. It does this by fetching the
    typedef and includes from the scheme-database.

    Example:
        A generated grid module::
            #include <dune/ufl/femscheme.hh>

            typedef FemScheme< DiffusionModel< Dune::Fem::AdaptiveLeafGridPart< Dune::YaspGrid< 2,
            Dune::EquidistantCoordinates< double, 2 > > >, 2 > , 2, fem > SchemeType;

            void registerDuneScheme ();

            BOOST_PYTHON_MODULE( dunescheme649df434446fad25fffcf9739ffb66d2 ) { registerDuneScheme(); }

    This would correspond to calling get("FemScheme", grid2d, 2), where grid2d is a python grid module.
    """
    return myGenerator.getModule(scheme, **parameters)

def get(scheme, space, grid, dimR, **parameters):
    """Call getModule() by passing in keyword arguments.
    """
    return getModule(scheme, space=space._module._typeName, gridpart=grid._module._typeName, dimRange=dimR, **parameters)

def create(scheme, space, grid, model, name, **parameters):
    """Get a Scheme.

    Call get() and create a C++ scheme class (see dune/fempy/dunescheme.hh).

    Args:
        scheme (string): the identifier for the scheme type to use
        space (Space): a space class generated from dune.fem.space
        grid (LeafGrid): a LeafGrid class generated from dune.fem.grid
        model (DiffusionModel) : a model class generated from dune.models.femufl
        name (string): the python name for the scheme
        parameters (kwargs): parameters used for fixing the scheme type

            * polorder=int: order of polynomials
            * solver=string: type of solver
    Returns:
        Scheme: the constructed scheme
    """
    module = get(scheme, space, grid, model.getDimRange(), **parameters)
    if hasattr(model, 'wrap'):
        scheme = module.Scheme(space, grid, model.wrap(), name)
    else:
        scheme = module.Scheme(space, grid, model, name)
    return scheme
