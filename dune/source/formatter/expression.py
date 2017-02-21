from __future__ import division, print_function, unicode_literals

from ...common.compatibility import isString

from ..builtin import *
from ..common import *
from ..expression import *
from ..operator import *


def entangleLists(lists, delimiter=''):
    left = lists[0]
    if len(lists) > 1:
        right = lists[1]
        n = len(left)
        return entangleLists([left[:n-1] + [left[n-1] + delimiter + right[0]] + right[1:]] + lists[2:], delimiter)
    else:
        return left


class FormatExpression:
    def __call__(self, expr):
        if isinstance(expr, Application):
            return self.application(expr)
        elif isinstance(expr, ConditionalExpression):
            return self.conditional(expr)
        elif isinstance(expr, ConstantExpression):
            return self.constant(expr)
        elif isinstance(expr, ConstructExpression):
            return self.construct(expr)
        elif isinstance(expr, DereferenceExpression):
            return self.dereference(expr)
        elif isinstance(expr, InitializerList):
            return self.initializerList(expr)
        elif isinstance(expr, LambdaExpression):
            return self.lambda_(expr)
        elif isinstance(expr, NullPtr):
            return self.nullptr(expr)
        elif isinstance(expr, Variable):
            return self.variable(expr)
        elif isinstance(expr, UnformattedExpression):
            return self.unformatted(expr)
        elif isString(expr):
            return self.unformatted(UnformattedExpression('auto', expr))
        else:
            raise Exception('Invalid type of expression: ' + str(type(expr)))

    def application(self, expr):
        args = [formatExpression(arg) for arg in expr.args]
        if isinstance(expr.function, Operator):
            if len(args) != expr.function.numArgs:
                raise Exception('The operator' + expr.function.name + ' takes ' + expr.function.numArgs + ' arguments (' + str(len(args)) + ' given).')
            if isinstance(expr.function, PrefixUnaryOperator):
                return entangleLists([[expr.function.name + '('], args[0], [')']])
            elif isinstance(expr.function, PostfixUnaryOperator):
                return entangleLists([['('], args[0], [')' + expr.function.name]])
            elif isinstance(expr.function, BinaryOperator):
                return entangleLists([['('], args[0], [') ' + expr.function.name + ' ('], args[1], [')']])
            elif isinstance(expr.function, BracketOperator):
                return entangleLists([['('], args[0], [')[ '], args[1], [' ]']])
            else:
                raise Exception('Unknown operator: ' + repr(expr.function))
        if isinstance(expr.function, BuiltInFunction):
            if expr.function.namespace is None:
                function = expr.function.name
            else:
                function = expr.function.namespace + '::' + expr.function.name
        elif isinstance(expr.function, (Function, Method)):
            function = expr.function.name
        else:
            function = expr.function
        if expr.args:
            return entangleLists([[function + '( '], entangleLists(args, ', '), [' )']])
        else:
            return [function + '()']

    def conditional(self, expr):
        return entangleLists([['('], self(expr.cond), [' ? '], self(expr.true), [' : '], self(expr.false), [')']])

    def constant(self, expr):
        return [expr.value]

    def construct(self, expr):
        if expr.args is not None:
            return entangleLists([[expr.cppType + '( '], entangleLists([self(arg) for arg in expr.args], ', '), [' )']])
        else:
            return [expr.cppType + '()']

    def dereference(self, expr):
        return entangleLists([['*'], self(expr.expr)])

    def initializerList(self, expr):
        return entangleLists([['{ '], entangleLists([self(arg) for arg in expr.args], ', '), [' }']])

    def lambda_(self, expr):
        capture = '' if not expr.capture else ' ' + ', '.join([c.name for c in expr.capture]) + ' '
        args = '' if expr.args is None else ' ' + ', '.join(expr.args) + ' '
        return ['[' + capture + '] (' + args + ') {', expr.code, '}']

    def nullptr(self, expr):
        return ['nullptr']

    def variable(self, expr):
        return [expr.name]

    def unformatted(self, expr):
        return [expr.value]


def formatExpression(expr):
    lines = FormatExpression()(expr)
    return lines
