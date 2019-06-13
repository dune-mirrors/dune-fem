from __future__ import division, print_function

from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.corealg.map_dag import map_expr_dags
from ufl.corealg.multifunction import MultiFunction
from ufl.differentiation import Grad


class DerivativeGatherer(MultiFunction):
    def __init__(self, arguments):
        MultiFunction.__init__(self)
        self.derivatives = {arg : set() for arg in arguments}

    def argument(self, expr):
        try:
            self.derivatives[expr].add(expr)
        except KeyError:
            pass

    def coefficient(self, expr):
        try:
            self.derivatives[expr].add(expr)
        except KeyError:
            pass

    def grad(self, expr):
        arg = expr
        while isinstance(arg, Grad):
            arg = arg.ufl_operands[0]
        try:
            self.derivatives[arg].add(expr)
        except KeyError:
            pass

    def expr(self, expr, *args):
        pass



def gatherDerivatives(form, arguments=None):
    if arguments is None:
        arguments = form.arguments()

    form = apply_derivatives(apply_algebra_lowering(form))

    gatherer = DerivativeGatherer(arguments)
    map_expr_dags(gatherer, [i.integrand() for i in form.integrals()])
    return [sorted(list(gatherer.derivatives[arg]), key=lambda d : len(d.ufl_shape)) for arg in arguments]
