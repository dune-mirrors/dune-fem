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
from . import discretefunction

myGenerator = generator.Generator("Scheme")

def solve( scheme, rhs=None, target=None, name=None, assemble=True ):
    if name == None:
        if hasattr(scheme, 'name'):
            name = scheme.name
        else:
            name == "no name"
    if target == None:
        if scheme.target == None:
            target = discretefunction.create(scheme._storage, scheme.space, name=name)
            target.interpolate( [0,]*scheme.dimRange )
        else:
            target = scheme.target
    if rhs == None:
        scheme._prepare()
    else:
        scheme._prepare(rhs)
    scheme._solve(target,assemble)
    return target

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
    module = myGenerator.getModule(scheme, **parameters)
    setattr(module.Scheme, "solve", solve)
    return myGenerator.getModule(scheme, **parameters)

def get(scheme, space, **parameters):
    """Call getModule() by passing in keyword arguments.
    """
    storage = parameters.get('storage', "Adaptive")
    try:
      nr = 65
      for s in space:
        discretefunction.get(storage, s._module, **parameters)
        parameters['space'+chr(nr)] = s._module._typeName
        nr += 1
      return getModule(scheme, **parameters)
    except:
      discretefunction.get(storage, space._module, **parameters)
      return getModule(scheme, space=space._module._typeName, **parameters)

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
        scheme = module.Scheme(grid, model.wrap(), name)
    else:
        scheme = module.Scheme(grid, model, name)
    return scheme
