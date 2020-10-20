from __future__ import absolute_import, division, print_function

from ufl.algorithms import expand_indices
from ufl.algorithms.analysis import extract_arguments
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.algorithms.transformer import Transformer
from ufl.constantvalue import IntValue, Zero
from ufl.differentiation import Grad
from ufl.restriction import Restricted

from .applyrestrictions import applyRestrictions

from .tensors import conditionalExprTensor, ExprTensor, keys

def sumTensorMaps(left, right):
    result = {key: l.copy() for key, l in left.items()}
    for key, r in right.items():
        l = left.get(key)
        result[key] = r.copy() if l is None else l + r
    return result


class MultiLinearExprSplitter(Transformer):
    def __init__(self, arguments):
        Transformer.__init__(self)
        self.arguments = tuple(arguments)
        self.empty = tuple(None for arg in arguments)

    def argument(self, expr):
        if len(expr.ufl_shape) > 0:
            raise Exception('Arguments should only occur in fully indexed expressions.')
        key = self._key(expr)
        if key != self.empty:
            return self._tensor(key, tuple())
        else:
            return self.terminal(expr)

    def indexed(self, expr):
        if len(expr.ufl_shape) > 0:
            raise Exception('Indexed expressions should always be scalar.')
        if len(expr.ufl_operands) != 2:
            raise Exception('Indexed expressions must have exactly two children.')
        operand = expr.ufl_operands[0]
        index = expr.ufl_operands[1]
        key = self._key(operand)
        if key != self.empty:
            return self._tensor(key, index.indices())
        else:
            return self.terminal(expr)

    def division(self, expr, left, right):
        if list(right.keys()) != [self.empty]:
            raise Exception('Only the left child of a division may access the linear arguments.')
        r = right[self.empty]
        return {key: l / r for key, l in left.items()}

    def negative_restricted(self, expr, arg):
        return {key: value.negative_restricted() for key, value in arg.items()}

    def positive_restricted(self, expr, arg):
        return {key: value.positive_restricted() for key, value in arg.items()}

    def product(self, expr, left, right):
        def oneOf(l, r):
            if l is None:
                return r
            if r is None:
                return l
            raise Exception('A linear argument may only appear in one child expression.')

        result = dict()
        for kl, el in left.items():
            for kr, er in right.items():
                key = tuple(oneOf(l, r) for l, r in zip(kl, kr))
                tensor = ExprTensor(self._shape(key))
                for i, il, ir in self._productIndices(kl, kr):
                    tensor[i] = el[il] * er[ir]
                result[key] = tensor
        return result

    def sum(self, expr, left, right):
        return sumTensorMaps(left, right)

    def conditional(self, expr):
        condition = expr.ufl_operands[0]
        trueCase = self.visit(expr.ufl_operands[1])
        falseCase = self.visit(expr.ufl_operands[2])

        result = dict()
        for key in set(trueCase.keys()) | set(falseCase.keys()):
            zero = ExprTensor(self._shape(key))
            result[key] = conditionalExprTensor(condition, trueCase.get(key, zero), falseCase.get(key, zero))
        return result

    def conj(self, expr, x):
        return x

    def max_value(self, expr, left, right):
        result = dict()
        if list(left.keys()) != [self.empty] or list(right.keys()) != [self.empty]:
            raise Exception('Linear arguments may not occur in maximum.')
        return { self.empty: expr }

    def min_value(self, expr, left, right):
        result = dict()
        if list(left.keys()) != [self.empty] or list(right.keys()) != [self.empty]:
            raise Exception('Linear arguments may not occur in minimum.')
        return { self.empty: expr }

    def terminal(self, expr):
        if len(expr.ufl_shape) > 0:
            raise Exception('Terminal expressions should only occur in fully indexed expressions.')
        tensor = ExprTensor(tuple())
        tensor[tuple()] = expr
        return {self.empty: tensor}

    atan = terminal
    atan_2 = terminal
    abs = terminal
    cos = terminal
    cosh = terminal
    exp = terminal
    ln = terminal
    power = terminal
    sin = terminal
    sinh = terminal
    sqrt = terminal
    tan = terminal
    tanh = terminal

    def variable(self, expr):
        return self.visit(expr.expression())

    def _getLeaf(self, expr):
        if isinstance(expr, Restricted):
            expr = expr.ufl_operands[0]
        while isinstance(expr, Grad):
            expr = expr.ufl_operands[0]
        return expr

    def _key(self, expr):
        leaf = self._getLeaf(expr)
        return tuple(expr if arg == leaf else None for arg in self.arguments)

    def _shape(self, key):
        shape = tuple()
        for arg in key:
            if arg is not None:
                shape += arg.ufl_shape
        return shape

    def _tensor(self, key, indices):
        tensor = ExprTensor(self._shape(key))
        tensor[indices] = IntValue(1)
        return {key: tensor}

    def _productIndices(self, keyl, keyr):
        if len(keyl) > 0:
            indices = self._productIndices(keyl[1:], keyr[1:])
            if keyl[0] is not None:
                key = keys(keyl[0].ufl_shape)
                return [(k+i, k+il, ir) for k in key for i, il, ir in indices]
            elif keyr[0] is not None:
                key = keys(keyr[0].ufl_shape)
                return [(k+i, il, k+ir) for k in key for i, il, ir in indices]
            return indices
        else:
            return [(tuple(), tuple(), tuple())]


def splitMultiLinearExpr(expr, arguments=None):
    if arguments is None:
        arguments = extract_arguments(expr)
    return MultiLinearExprSplitter(arguments).visit(expr)



def splitForm(form, arguments=None):
    if arguments is None:
        arguments = form.arguments()

    form = applyRestrictions(form)
    form = expand_indices(apply_derivatives(apply_algebra_lowering(form)))
    form = applyRestrictions(form)

    integrals = {}
    for integral in form.integrals():
        right = splitMultiLinearExpr(integral.integrand(), arguments)
        left = integrals.get(integral.integral_type())
        if left is not None:
            right = sumTensorMaps(left, right)
        integrals[integral.integral_type()] = right

    return integrals
