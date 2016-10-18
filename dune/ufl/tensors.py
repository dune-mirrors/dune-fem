from __future__ import division, print_function

from ufl import as_tensor
from ufl.constantvalue import Zero
from ufl.core.multiindex import FixedIndex, MultiIndex


# ExprTensor
# ----------

class ExprTensor:
    def __init__(self, shape, data=None):
        self.shape = shape
        self.data = self._zero(shape) if data is None else data

    def __repr__(self):
        return repr(self.data)

    def __add__(self, other):
        if not isinstance(other, ExprTensor):
            raise Exception('Cannot add ' + type(other) + ' to ' + type(self) + '.')
        if other.shape != self.shape:
            raise Exception('Cannot add tensors of different shape.')
        return ExprTensor(self.shape, self._add(self.shape, self.data, other.data))

    def __truediv__(self, other):
        if isinstance(other, ExprTensor):
            raise Exception('Cannot divide by tensors tensors.')
        return ExprTensor(self.shape, self._div(self.shape, self.data, other))

    def __mul__(self, other):
        if isinstance(other, ExprTensor):
            raise Exception('Cannot multiply tensors.')
        return ExprTensor(self.shape, self._mul(self.shape, self.data, other))

    def __getitem__(self, key):
        if isinstance(key, MultiIndex):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(shape)) + '.')
            if not all(isinstance(idx, FixedIndex) for idx in key):
                raise Exception('Expect only FixedIndex entries in MultiIndex.')
            data = self.data
            key = key.indices()
            while len(key) > 0:
                data = data[int(key[0])]
                key = key[1:]
            return data
        else:
            raise Exception('Expect MultiIndex to access component')

    def __setitem__(self, key, value):
        if isinstance(key, MultiIndex):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(shape)) + '.')
            if not all(isinstance(idx, FixedIndex) for idx in key):
                raise Exception('Expect only FixedIndex entries in MultiIndex.')
            key = key.indices()
            if len(key) > 0:
                data = self.data
                while len(key) > 1:
                    data = data[int(key[0])]
                    key = key[1:]
                data[int(key[0])] = value
            else:
                self.data = value
        else:
            raise Exception('Expect MultiIndex to access component')

    def keys(self):
        return [MultiIndex(idx) for idx in self._keys(self.shape)]

    def as_ufl(self):
        return self._as_ufl(self.shape, self.data)

    def _add(self, shape, left, right):
        if len(shape) > 0:
            return [self._add(shape[1:], left[i], right[i]) for i in range(0, shape[0])]
        else:
            return left + right

    def _as_ufl(self, shape, data):
        if len(shape) > 0:
            return as_tensor([self._as_ufl(shape[1:], data[i]) for i in range(0, shape[0])])
        else:
            return data

    def _keys(self, shape):
        if len(shape) > 0:
            return [(FixedIndex(i),) + t for i in range(0, shape[0]) for t in self._keys(shape[1:])]
        else:
            return tuple()

    def _div(self, shape, tensor, value):
        if len(shape) > 0:
            return [self._div(shape[1:], tensor[i], value) for i in range(0, shape[0])]
        else:
            return tensor / value

    def _mul(self, shape, tensor, value):
        if len(shape) > 0:
            return [self._mul(shape[1:], tensor[i], value) for i in range(0, shape[0])]
        else:
            return tensor * value

    def _zero(self, shape):
        if len(shape) > 0:
            return [self._zero(shape[1:]) for i in range(0, shape[0])]
        else:
            return Zero()
