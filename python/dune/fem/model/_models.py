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
    # the following needs to be done for the linearized problem:
    # if isinstance(equation, ufl.equation.Equation):
    #     lhs = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(equation.lhs)))
    #     if lhs == ufl.adjoint(lhs):
    #         setattr(Model, 'symmetric', 'true')
    #    else:
    #        setattr(Model, 'symmetric', 'false')

    return Model(coefficients=coefficients)

def integrands(view, form, *args, **kwargs):
    import ufl
    import dune.ufl
    import dune.models.integrands as integrands

    tempVars   = kwargs.pop('tempVars', True)
    virtualize = kwargs.pop('virtualize',True)
    modelPatch = kwargs.pop('modelPatch',None)
    includes   = kwargs.pop('includes',None)

    Integrands = integrands.load(view, form, *args,
                     tempVars=tempVars,virtualize=virtualize,
                     modelPatch=modelPatch,includes=includes)

    return Integrands(*args, **kwargs)
