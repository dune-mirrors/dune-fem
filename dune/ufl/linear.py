from __future__ import absolute_import, division, print_function

from ufl.algorithms.analysis import extract_arguments
from ufl.algorithms.transformer import Transformer
from ufl.constantvalue import IntValue
from ufl.differentiation import Grad

from .tensors import ExprTensor, keys

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
            tensor = ExprTensor(self._shape(key))
            tensor[tuple()] = IntValue(1)
            return {key: tensor}
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
            tensor = ExprTensor(operand.ufl_shape)
            tensor[index.indices()] = IntValue(1)
            return {key: tensor}
        else:
            return self.terminal(expr)

    def division(self, expr, left, right):
        if list(right.keys()) != [self.empty]:
            raise Exception('Only the left child of a division may access the linear arguments.')
        r = right[self.empty]
        return {key: l / r for key, l in left.items()}

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
        for key, r in right.items():
            l = left.get(key)
            left[key] = r if l is None else l + r
        return left

    def terminal(self, expr):
        if len(expr.ufl_shape) > 0:
            raise Exception('Terminal expressions should only occur in fully indexed expressions.')
        tensor = ExprTensor(tuple())
        tensor[tuple()] = expr
        return {self.empty: tensor}

    atan = terminal
    atan_2 = terminal
    cos = terminal
    sin = terminal
    power = terminal
    tan = terminal

    def variable(self, expr):
        return self.visit(expr.expression())

    def _getLeaf(self, expr):
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
