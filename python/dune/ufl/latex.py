from __future__ import absolute_import, print_function

from ufl.core.multiindex import FixedIndex, MultiIndex

from ufl.corealg.map_dag import map_expr_dag
from ufl.algorithms.transformer import Transformer

class Expr2Latex(Transformer):
    def __init__(self, arguments=None, coefficients=None):
        Transformer.__init__(self)
        self.arguments = arguments
        self.coefficients = coefficients

    def argument(self, expr):
        if self.arguments is None:
            return "\\varphi^{" + str(expr.number()) + "}"
        else:
            return arguments[expr]

    def atan(self, expr, arg):
        return "\\atan\\left(" + arg + "\\right)"

    def atan_2(self, expr, left, right):
        return "\\atan\\left(\\frac{" + left + "}{" + right + "}\\right)"

    def coefficient(self, expr):
        if self.coefficients is None:
            return "w^{" + str(expr.count()) + "}"
        else:
            return coefficients[expr]

    def cos(self, expr, arg):
        return "\\cos\\left(" + arg + "\\right)"

    def division(self, expr, left, right):
        return "\\frac{" + left + "}{" + right + "}"

    def dot(self, expr, left, right):
        return left + "\\cdot" + right

    def exp(self, expr, arg):
        return "\\exp\\left(" + arg + "\\right)"

    def float_value(self, expr):
        value = expr.value();
        if value < 0:
            return "(" + str(value) + ")"
        else:
            return str(value)

    def grad(self, expr, arg):
        return "\\nabla " + arg

    def indexed(self, expr, arg, index):
        return arg + "_" + index

    def inner(self, expr, left, right):
        return left + "\\cdot" + right

    int_value = float_value

    def list_tensor(self, expr):
        shape = expr.ufl_shape
        if len(shape) == 1:
            lines = [self.visit(op) for op in expr.ufl_operands]
        else:
            raise Exception("Output not implemented for tensor of rank " + str(len(shape)) + ".")
        return "\\begin{pmatrix}" + "\\\\ \n".join(lines) + "\\end{pmatrix}"

    def multi_index(self, expr):
        return "{" + ",".join([str(i) for i in expr]) + "}"

    def negative_restricted(self, expr, arg):
        return "\\left(" + arg + "\\right)^-"

    def positive_restricted(self, expr, arg):
        return "\\left(" + arg + "\\right)^+"

    def product(self, expr, left, right):
        return left + "\," + right

    def power(self, expr, base, exponent):
        return base + "^{" + exponent + "}"

    def sin(self, expr, arg):
        return "\\sin\\left(" + arg + "\\right)"

    def spatial_coordinate(self, expr):
        return "x"

    def facet_normal(self, expr):
        return "n"

    def sum(self, expr, left, right):
        return left + " + " + right

    def tan(self, expr, arg):
        return "\\tan\\left(" + arg + "\\right)"

    def zero(self, expr):
        return "0"


def expr2latex(expr):
    return Expr2Latex().visit(expr)


def form2latex(form):
    result = ""
    for integral in form.integrals():
        if integral.integral_type() == 'cell':
            result += "\\int_{\\Omega}" + expr2latex(integral.integrand()) + "\,dx"
        elif integral.integral_type() == 'exterior_facet':
            result += "\\int_{\\partial\\Omega}" + expr2latex + "\,dx"
    return result


def equation2latex(equation):
    try:
        return form2latex(equation.lhs) + "=" + form2latex(equation.rhs)
    except AttributeError:
        return form2latex(equation) + "= 0"
