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


def formatExpression(expr):
    if isinstance(expr, Application):
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
    elif isinstance(expr, ConditionalExpression):
        return entangleLists([['('], formatExpression(expr.cond), [' ? '], formatExpression(expr.true), [' : '], formatExpression(expr.false), [')']])
    elif isinstance(expr, ConstantExpression):
        return [expr.value]
    elif isinstance(expr, ConstructExpression):
        if expr.args is not None:
            return entangleLists([[expr.cppType + '( '], entangleLists([formatExpression(arg) for arg in expr.args], ', '), [' )']])
        else:
            return [expr.cppType + '()']
    elif isinstance(expr, DereferenceExpression):
        return entangleLists([['*'], formatExpression(expr.expr)])
    elif isinstance(expr, InitializerList):
        return entangleLists([['{ '], entangleLists([formatExpression(arg) for arg in expr.args], ', '), [' }']])
    elif isinstance(expr, LambdaExpression):
        capture = '' if not expr.capture else ' ' + ', '.join([c.name for c in expr.capture]) + ' '
        args = '' if expr.args is None else ' ' + ', '.join(expr.args) + ' '
        return ['[' + capture + '] (' + args + ') {', expr.code, '}']
    elif isinstance(expr, NullPtr):
        return ['nullptr']
    elif isinstance(expr, Variable):
        return [expr.name]
    elif isinstance(expr, UnformattedExpression):
        return [expr.value]
    elif isString(expr):
        return [expr.strip()]
    else:
        raise Exception('Invalid type of expression: ' + str(type(expr)))
