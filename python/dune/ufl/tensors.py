from __future__ import division, print_function

from ufl import as_tensor
from ufl.conditional import Conditional
from ufl.constantvalue import Zero
from ufl.core.multiindex import MultiIndex
from ufl.restriction import NegativeRestricted, PositiveRestricted


def keys(shape):
    if len(shape) > 0:
        return [(i,) + t for i in range(0, shape[0]) for t in keys(shape[1:])]
    else:
        return [tuple()]


def apply(op, shape, *args):
    if False: #len(shape) > 2:
        return [apply(op, shape[1:], *(arg[i,:,:] for arg in args)) for i in range(0, shape[0])]
    elif False: #len(shape) > 1:
        return [apply(op, shape[1:], *(arg[i,:] for arg in args)) for i in range(0, shape[0])]
    elif len(shape) > 0:
        return [apply(op, shape[1:], *(arg[i] for arg in args)) for i in range(0, shape[0])]
    else:
        return op(*args)


def reformat(fmt, shape, data):
    if len(shape) > 0:
        return fmt([reformat(fmt, shape[1:], data[i]) for i in range(0, shape[0])])
    else:
        return data


def fill(shape, value):
    if len(shape) > 0:
        return [fill(shape[1:], value) for i in range(0, shape[0])]
    else:
        return value


def getItem(data, key):
    while len(key) > 0:
        data = data[int(key[0])]
        key = key[1:]
    return data


def setItem(data, key, value):
    result = data
    if len(key) > 0:
        while len(key) > 1:
            data = data[int(key[0])]
            key = key[1:]
        data[int(key[0])] = value
        return result
    else:
        return value



# ExprTensor
# ----------

class ExprTensor:
    def __init__(self, shape, data=None):
        self.shape = shape
        self.data = fill(shape, Zero()) if data is None else data

    def __repr__(self):
        return repr(self.data)

    def __add__(self, other):
        if not isinstance(other, ExprTensor):
            raise Exception('Cannot add ' + type(other) + ' to ' + type(self) + '.')
        if other.shape != self.shape:
            raise Exception('Cannot add tensors of different shape.')
        return ExprTensor(self.shape, apply(lambda u, v: u + v, self.shape, self.data, other.data))

    def __truediv__(self, other):
        if isinstance(other, ExprTensor):
            if len(other.shape) > 0:
                raise Exception('Cannot divide by tensors tensors.')
            else:
                other = other.data
        return ExprTensor(self.shape, apply(lambda v: v / other, self.shape, self.data))

    def __mul__(self, other):
        if isinstance(other, ExprTensor):
            raise Exception('Cannot multiply tensors.')
        return ExprTensor(self.shape, apply(lambda v: v * other, self.shape, self.data))

    def __neg__(self):
        return ExprTensor(self.shape, apply(lambda v: -v, self.shape, self.data))

    def __getitem__(self, key):
        if isinstance(key, MultiIndex):
            key = key.indices()
        if isinstance(key, tuple):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(self.shape)) + ' (got ' + str(len(key)) + ').')
            return getItem(self.data, key)
        else:
            raise Exception('Expect tuple or MultiIndex to access component')

    def __setitem__(self, key, value):
        if isinstance(key, MultiIndex):
            key = key.indices()
        if isinstance(key, tuple):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(self.shape)) + ' (got ' + str(len(key)) + ').')
            self.data = setItem(self.data, key, value)
        else:
            raise Exception('Expect tuple or MultiIndex to access component')

    def copy(self):
        return ExprTensor(self.shape, apply(lambda v: v, self.shape, self.data))

    def negative_restricted(self):
        return ExprTensor(self.shape, apply(lambda v: NegativeRestricted(v), self.shape, self.data))

    def positive_restricted(self):
        return ExprTensor(self.shape, apply(lambda v: PositiveRestricted(v), self.shape, self.data))

    def for_each(self, op):
        return ExprTensor(self.shape, apply(lambda v: op(v), self.shape, self.data))

    def keys(self):
        return keys(self.shape)

    def as_ufl(self):
        return reformat(as_tensor, self.shape, self.data)

    def is_zero(self):
        return self.as_ufl() == ExprTensor(self.shape).as_ufl()


def conditionalExprTensor(condition, trueTensor, falseTensor):
    if not isinstance(trueTensor, ExprTensor) or not isinstance(falseTensor, ExprTensor):
        raise Exception('conditionalExprTensor works on ExprTensors only.')
    if trueTensor.shape != falseTensor.shape:
        raise Exception('Cannot construct conditional for tensors of different shape.')
    shape = trueTensor.shape
    def f(u, v):
        if isinstance(u, Zero) and isinstance(v, Zero):
            return u
        else:
            return Conditional(condition, u, v)
    return ExprTensor(shape, apply(f, shape, trueTensor.data, falseTensor.data))
