from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import inspect

from dune.generator.generator import SimpleGenerator

from ._functions import *

import dune.ufl

def integrate(grid,expression,order):
    try:
        return expression.integrate()
    except AttributeError:
        return uflFunction(grid,"tmp",order,expression).integrate()
    # return dune.ufl.expression2GF(grid,expression,order).integrate()
# perhaps a general
#    def assemble(grid/space, expression, order):
# would be better. If expression is a function this would return a
# fieldvector with the integral of this function, if expression is a
# functional a discreteFunctional and for a bilinear form a matrix is
# returned. This could also include a "piecewise" option that (at least in
# the case of a function) returns the integral on each element.
