from __future__ import division, print_function, unicode_literals

from .model import Integrands
from .ufl import compileUFL
from .load import load

def create(grid, integrands, renumbering=None, tempVars=True):
    return load(grid, integrands, renumbering=renumbering, tempVars=tempVars).Integrands()
