"""Functions for creating python modules and C++ classes for schemes.
"""

from __future__ import print_function
import importlib
import hashlib
import os.path
import re
# from termcolor import colored

from dune.deprecate import deprecated
from ._schemes import *

from dune.fem.deprecated import deprecated
from dune.generator.generator import SimpleGenerator

def solve( scheme, target, rhs=None, *, rightHandSide=None ):
    if rhs is not None:
        deprecated("""
The `solve` method with argument `rhs` is deprecated, used named argument ``rightHandSide`` instead.
Note that the behavior has changed if dirichlet boundary constraints are present.
When using `rhs` argument the result on the boundary is `target=g`
while the new behavior leads to `target=rightHandSide+g`.

See changelog entry in tutorial for more details.
""")
        rightHandSide = rhs.copy()
        try:
            scheme.setConstraints(0,rightHandSide) # this is what was implemented originally in the scheme
        except:
            pass # no constraints available

    if rightHandSide is None:
        rightHandSide = scheme.space.zero

    if hasattr(scheme,"preconditioning") and scheme.preconditioning is not None:
        try:
            # if this call fails then because the SFINAE in fempy/py/scheme.h has failed
            return scheme._solve(solution=target,
                                 rightHandSide=rightHandSide,
                                 preconditioning=scheme.preconditioning)
        except:
            raise Exception("scheme.solve with Python based preconditioning failed!")
    else:
        return scheme._solve(solution=target, rightHandSide=rightHandSide)

_defaultGenerator = SimpleGenerator("Scheme", "Dune::FemPy")

def addAttr(module, cls):
    setattr(cls, "solve", solve)

fileBase = "femscheme"

def module(includes, typeName, *args, backend=None,
           generator=_defaultGenerator,
           baseClasses=None):
    from dune.fem.space import addBackend
    includes = includes + ["dune/fempy/py/scheme.hh"]
    moduleName = fileBase + "_" + hashlib.md5(typeName.encode('utf-8')).hexdigest()
    module = generator.load(includes, typeName, moduleName, *args, dynamicAttr=True, baseClasses=baseClasses)
    addAttr(module, module.Scheme)
    JacobianOperator = getattr(module.Scheme,"JacobianOperator",None)
    if JacobianOperator is not None and hasattr(JacobianOperator,"_backend") and backend is not None:
        addBackend(JacobianOperator,backend)
    return module
