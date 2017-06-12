from __future__ import absolute_import, division, print_function, unicode_literals
__metaclass__ = type

import hashlib
import importlib
import inspect

from dune.generator.generator import SimpleGenerator

from ._functions import *

def integrate(grid,expression,order):
    shape = expression.ufl_shape
    assert len(shape) == 0 or len(shape) == 1 , "can only integrate scalar or vector valued expression"
    if len(shape) == 0:
        gf = localFunction(grid, "tmp", order, lambda e,x: [expression(e(x))] )
        return gf.integrate()[0]
    if len(shape) == 1:
        gf = localFunction(grid, "tmp", order, lambda e,x: [ expression[i](e(x)) for i in range(shape[0]) ] )
        return gf.integrate()
# perhaps a general
#    def assemble(grid/space, expression, order):
# would be better. If expression is a function this would return a
# fieldvector with the integral of this function, if expression is a
# functional a discreteFunctional and for a bilinear form a matrix is
# returned. This could also include a "piecewise" option that (at least in
# the case of a function) returns the integral on each element.
