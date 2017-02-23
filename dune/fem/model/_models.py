from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

# def elliptic(view, equation, dirichlet = {}, exact = None, tempVars=True, coefficients={}, header=False):
def elliptic(view, equation, *args, **kwargs):
    import ufl
    import dune.ufl
    import dune.models.elliptic as elliptic
    Model = elliptic.importModel(view, equation, tempVars=False, *args, **kwargs).Model
    coefficients = kwargs.pop("coefficients", {})
    if isinstance(equation, ufl.equation.Equation):
        form = equation.lhs - equation.rhs
        uflCoeff = set(form.coefficients())
        # for bndId in dirichlet:
        #     for expr in dirichlet[bndId]:
        #         _, c = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
        #         uflCoeff |= set(c)
        fullCoeffs = {c:c.gf for c in uflCoeff if isinstance(c,dune.ufl.GridCoefficient)}
        lhs = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(equation.lhs)))
        if lhs == ufl.adjoint(lhs):
            setattr(Model, 'symmetric', 'true')
        else:
            setattr(Model, 'symmetric', 'false')
    else:
        fullCoeffs = {}
    print(fullCoeffs, coefficients)
    fullCoeffs.update(coefficients)
    return Model( coefficients=fullCoeffs )


def integrands(view, equation, tempVars=True, coefficients={}):
    import ufl
    import dune.ufl
    import dune.models.integrands as integrands
    return integrands.create(view, equation, tempVars=tempVars)
