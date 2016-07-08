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
            name == "default"
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
            #include <dune/fempy/py/scheme.hh>

            #include <dune/alugrid/grid.hh>
            #include <dune/alugrid/dgf.hh>
            #include <dune/fem/gridpart/adaptiveleafgridpart.hh>
            #include <dune/fem/space/lagrange.hh>
            #include <dune/fem/schemes/femscheme.hh>

            typedef FemScheme< ... > DuneType;

            PYBIND11_PLUGIN( scheme_325fb4c39707a927cfd2ab14f477d3fc )
            {
              pybind11::module module( "scheme_325fb4c39707a927cfd2ab14f477d3fc" );
              Dune::FemPy::registerScheme< DuneType >( module );
              return module.ptr();
            }

    This would correspond to calling get("FemScheme", space), where space is a python space module.
    """
    module = myGenerator.getModule(scheme, **parameters)
    setattr(module.Scheme,"solve", solve)
    setattr(module.Scheme, "_storage", module._selector.parameters["storage"])
    setattr(module.Scheme, "target", None)
    return module

def get(scheme, space, **parameters):
    """Call getModule() by passing in keyword arguments.
    """
    storage = parameters.pop('storage', "Adaptive")
    try:
      nr = 65    # look for spaceA,spaceB,...
      extra_includes=""
      for s in space:
        dfmodule = discretefunction.get(storage, s._module, **parameters)
        storage  = dfmodule.DiscreteFunction._storage
        parameters['space'+chr(nr)] = s._module._typeName
        extra_includes += s._module._includes + dfmodule._includes
        nr += 1
    except:
      dfmodule = discretefunction.get(storage, space._module, **parameters)
      storage  = dfmodule.DiscreteFunction._storage
      parameters['space'] = space._module._typeName
      extra_includes=space._module._includes + dfmodule._includes

    module = getModule(scheme,
                       storage=storage,
                       extra_includes=extra_includes,
                       **parameters)
    return module

def create(scheme, space_or_target, model, name, *param, **parameters):
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
    try:
        space = space_or_target.space
        target = space_or_target
    except:
        space = space_or_target
        target = None

    module = get(scheme, space, **parameters)

    class ExtendedScheme(module.Scheme):
        def __init__(self,space,model,name,target=None):
            module.Scheme.__init__(self,space,model,name, *param)
            self.target = target

    ret = ExtendedScheme(space, model, name, target)

    return ret
