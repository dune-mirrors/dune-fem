from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

def elliptic(view, equation, *args, **kwargs):
    import ufl
    import dune.ufl
    import dune.models.elliptic as elliptic

    coefficients = kwargs.pop('coefficients', dict())

    Model = elliptic.load(view, equation, *args, **kwargs).Model
    if isinstance(equation, ufl.equation.Equation):
        lhs = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(equation.lhs)))
        if lhs == ufl.adjoint(lhs):
            setattr(Model, 'symmetric', 'true')
        else:
            setattr(Model, 'symmetric', 'false')

    return Model(coefficients=coefficients)


def integrands(view, equation, tempVars=True, coefficients={}):
    import ufl
    import dune.ufl
    import dune.models.integrands as integrands
    return integrands.create(view, equation, tempVars=tempVars)
