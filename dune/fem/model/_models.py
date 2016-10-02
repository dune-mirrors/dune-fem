from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration

def elliptic(grid, equation, dirichlet = {}, exact = None, tempVars=True, coefficients={}):
    import ufl
    import dune.models.elliptic as elliptic
    Model = elliptic.importModel(grid, equation, dirichlet, exact, tempVars).Model
    if isinstance(equation, ufl.equation.Equation):
        form = equation.lhs - equation.rhs
        uflCoeff = set(form.coefficients())
        for bndId in dirichlet:
            for expr in dirichlet[bndId]:
                _, c = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
                uflCoeff |= set(c)
        fullCoeffs = {c:c.gf for c in uflCoeff if isinstance(c,GridCoefficient)}
    else:
        fullCoeffs = {}
    fullCoeffs.update(coefficients)
    return Model( coefficients=fullCoeffs )
