from __future__ import division, print_function, unicode_literals

from types import GeneratorType

from ..common.utility import isInteger

from .common import Block, Statement
from .operator import BinaryOperator, BracketOperator, PrefixUnaryOperator, PostfixUnaryOperator


class Expression(Statement):
    def __init__(self, cppType=None):
        Statement.__init__(self)
        self.cppType = cppType

    def __add__(self, other):
        return Application(BinaryOperator('+'), args=(self, makeExpression(other)))

    def __sub__(self, other):
        return Application(BinaryOperator('-'), args=(self, makeExpression(other)))

    def __mul__(self, other):
        return Application(BinaryOperator('*'), args=(self, makeExpression(other)))

    def __truediv__(self, other):
        return Application(BinaryOperator('/'), args=(self, makeExpression(other)))

    def __mod__(self, other):
        return Application(BinaryOperator('%'), args=(self, makeExpression(other)))

    def __iadd__(self, other):
        return Application(BinaryOperator('+='), args=(self, makeExpression(other)))

    def __isub__(self, other):
        return Application(BinaryOperator('-='), args=(self, makeExpression(other)))

    def __imul__(self, other):
        return Application(BinaryOperator('*='), args=(self, makeExpression(other)))

    def __imod__(self, other):
        return Application(BinaryOperator('%='), args=(self, makeExpression(other)))

    def __itruediv__(self, other):
        return Application(BinaryOperator('/='), args=(self, makeExpression(other)))

    def __lt__(self, other):
        return Application(BinaryOperator('<'), args=(self, makeExpression(other)))

    def __le__(self, other):
        return Application(BinaryOperator('<='), args=(self, makeExpression(other)))

    def __eq__(self, other):
        return Application(BinaryOperator('=='), args=(self, makeExpression(other)))

    def __ne__(self, other):
        return Application(BinaryOperator('!='), args=(self, makeExpression(other)))

    def __ge__(self, other):
        return Application(BinaryOperator('>='), args=(self, makeExpression(other)))

    def __gt__(self, other):
        return Application(BinaryOperator('>'), args=(self, makeExpression(other)))

    def __neg__(self):
        return Application(PrefixUnaryOperator('-'), args=(self,))

    def __pos__(self):
        return Application(PrefixUnaryOperator('+'), args=(self,))

    def __getitem__(self, index):
        if self.cppType is not None and self.cppType.startswith('std::tuple'):
            from .builtin import get
            return Application(get(index), args=(self,))
        else:
            if isinstance(index, tuple):
                result = self
                for i in index:
                    result = Application(BracketOperator(), args=(result, makeExpression(i)))
                return result
            else:
                return Application(BracketOperator(), args=(self, makeExpression(index)))


class Application(Expression):
    def __init__(self, function, args=None):
        Expression.__init__(self)
        self.function = function
        self.args = tuple(args)

    def __hash__(self):
        return hash((self.cppType, self.function, self.args))


class ConditionalExpression(Expression):
    def __init__(self, cppType, cond, true, false):
        Expression.__init__(self, cppType)
        self.cond = cond
        self.true = true
        self.false = false

    def __hash__(self):
        return hash((self.cppType, self.cond, self.true, self.false))


class ConstantExpression(Expression):
    def __init__(self, cppType, value):
        Expression.__init__(self, cppType)
        self.value = value

    def __hash__(self):
        return hash((self.cppType, self.value))


class ConstructExpression(Expression):
    def __init__(self, cppType, args=None, brace=False):
        Expression.__init__(self, cppType)
        self.args = None if args is None else [makeExpression(arg) for arg in args]
        self.brace = brace

    def open(self):
        return '{' if self.brace else '('

    def close(self):
        return '}' if self.brace else ')'

    def __hash__(self):
        return hash((self.cppType, self.args, self.brace))


class DereferenceExpression(Expression):
    def __init__(self, expr):
        Expression.__init__(self, 'auto')
        self.expr = expr

    def __hash__(self):
        return hash('operator*', self.expr)


class InitializerList(Expression):
    def __init__(self, *args):
        Expression.__init__(self)
        self.args = tuple(args)


class LambdaExpression(Expression):
    def __init__(self, args=None, capture=None, code=None):
        Expression.__init__(self, None)
        self.args = None if args is None else tuple(args)
        if capture is None:
            self.capture = None
        else:
            self.capture = sorted(capture, key=lambda x: x.name)
        if code is None:
            self.code = None
        elif isinstance(code, Block):
            self.code = tuple(block.content)
        elif isinstance(code, (list, set, tuple, GeneratorType)):
            self.code = tuple(code)
        else:
            self.code = (code,)

    def __call__(self, *args):
        if (len(args) != len(self.args)):
            raise Exception('Wrong number of Arguments: ' + str(len(args)) + ' (should be ' + str(len(self.args)) + ').')
        return Application(self, args=[makeExpression(arg) for arg in args])

    def __hash__(self):
        return hash((self.cppType, self.args, self.capture, self.code))


class NullPtr(Expression):
    def __init__(self):
        Expression.__init__(self, "std::nullptr_t")

    def __hash__(self):
        return hash(self.cppType)


class This(Expression):
    def __init__(self):
        Expression.__init__(self, "auto")

    name = 'this'

    def __hash__(self):
        return hash(("auto", "this"))


class Variable(Expression):
    def __init__(self, cppType, name):
        Expression.__init__(self, cppType)
        self.name = name

    def __hash__(self):
        return hash((self.cppType, self.name))


class UnformattedExpression(Expression):
    def __init__(self, cppType, value, uses=None):
        Expression.__init__(self, cppType)
        self.value = value.strip()
        self.uses = uses

    def __hash__(self):
        return hash((self.cppType, self.value))


def makeExpression(expr):
    if isinstance(expr, Expression):
        return expr
    elif isinstance(expr, bool):
        return ConstantExpression('bool', 'true' if expr else 'false')
    elif isInteger(expr):
        if expr < 0:
            return -ConstantExpression('long', str(-expr))
        else:
            return ConstantExpression('long', str(expr))
    elif isinstance(expr, float):
        s = str(abs(expr))
        if ("." not in s) and ("e" not in s):
            s += ".0"
        e = ConstantExpression('double', s)
        return -e if expr < 0 else e
    else:
        # print("Warning: Making expression from " + str(expr) + "(type: " + str(type(expr)) + ").")
        return expr


def assign(left, right):
    return Application(BinaryOperator('='), args=(left, makeExpression(right)))


def construct(cppType, *args, **kwargs):
    return ConstructExpression(cppType, args if args else None, **kwargs)


def dereference(expr):
    return DereferenceExpression(expr)


def lambda_(args=None, capture=None, code=None):
    return LambdaExpression(args=args, capture=capture, code=code)


nullptr = NullPtr()
this = This()
