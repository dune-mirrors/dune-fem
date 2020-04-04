from __future__ import division, print_function, unicode_literals

from ...common.utility import isString

from ..builtin import *
from ..common import *
from ..expression import *
from ..function import *
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
    _priority_terminal = 0
    _priority_postfix = 1
    _priority_prefix = 2
    _priority_multiplication = 3
    _priority_addition = 4
    _priority_relational = 5
    _priority_equality = 6
    _priority_assignment = 7
    _priority_none = 8

    _priority_binary = {
        '+':  _priority_addition,
        '-':  _priority_addition,
        '*':  _priority_multiplication,
        '/':  _priority_multiplication,
        '%':  _priority_multiplication,
        '=':  _priority_assignment,
        '+=': _priority_assignment,
        '-=': _priority_assignment,
        '*=': _priority_assignment,
        '%=': _priority_assignment,
        '/=': _priority_assignment,
        '<':  _priority_relational,
        '<=': _priority_relational,
        '==': _priority_equality,
        '!=': _priority_equality,
        '>=': _priority_relational,
        '>':  _priority_relational
        }

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
        elif isinstance(expr, This):
            return self.this(expr)
        elif isinstance(expr, Variable):
            return self.variable(expr)
        elif isinstance(expr, UnformattedExpression):
            return self.unformatted(expr)
        elif isString(expr):
            return self.unformatted(UnformattedExpression('auto', expr))
        else:
            raise Exception('Invalid type of expression: ' + str(type(expr)))

    def formatArg(self, priority, arg):
        arg = self(arg)
        return arg[0] if priority > arg[1] else entangleLists([['('], arg[0], [')']])

    def application(self, expr):
        if isinstance(expr.function, Operator):
            return self.operatorApplication(expr)

        priority = FormatExpression._priority_postfix
        if isinstance(expr.function, BuiltInFunction):
            if expr.function.namespace is None:
                function = [expr.function.name]
            else:
                function = [expr.function.namespace + '::' + expr.function.name]
        elif isinstance(expr.function, (Function, Method)):
            function = [expr.function.name]
        elif isinstance(expr.function, LambdaExpression):
            function = self.formatArg(priority, expr.function)
        else:
            function = [expr.function]

        if expr.args:
            args = [self.formatArg(FormatExpression._priority_none, arg) for arg in expr.args]
            return entangleLists([function, ['( '], entangleLists(args, ', '), [' )']]), priority
        else:
            return entangleLists([function, ['()']]), priority

    def conditional(self, expr):
        priority = FormatExpression._priority_assignment
        cond = self.formatArg(priority, expr.cond)
        true = self.formatArg(priority, expr.true)
        false = self.formatArg(priority, expr.false)
        return entangleLists([cond, [' ? '], true, [' : '], false]), priority

    def constant(self, expr):
        return [expr.value], FormatExpression._priority_terminal

    def construct(self, expr):
        priority = FormatExpression._priority_postfix
        if expr.args is not None:
            args = [self.formatArg(FormatExpression._priority_none, arg) for arg in expr.args]
            return entangleLists([[expr.cppType + expr.open() + ' '], entangleLists(args, ', '), [' ' + expr.close()]]), priority
        else:
            return [expr.cppType + expr.open() + expr.close()], priority

    def dereference(self, expr):
        priority = FormatExpression._priority_prefix
        arg = self.formatArg(priority, expr.expr)
        return entangleLists([['*'], arg]), priority

    def initializerList(self, expr):
        args = [self.formatArg(FormatExpression._priority_none, arg) for arg in expr.args]
        return entangleLists([['{ '], entangleLists(args, ', '), [' }']]), FormatExpression._priority_terminal

    def lambda_(self, expr):
        capture = '' if not expr.capture else ' ' + ', '.join([c.name for c in expr.capture]) + ' '
        args = '' if expr.args is None else ' ' + ', '.join(expr.args) + ' '
        return ['[' + capture + '] (' + args + ') {', expr.code, '}'], FormatExpression._priority_terminal

    def nullptr(self, expr):
        return ['nullptr'], FormatExpression._priority_terminal

    def operatorApplication(self, expr):
        if len(expr.args) != expr.function.numArgs:
            raise Exception('The operator' + expr.function.name + ' takes ' + expr.function.numArgs + ' arguments (' + str(len(expr.args)) + ' given).')
        if isinstance(expr.function, PrefixUnaryOperator):
            return entangleLists([[expr.function.name], self.formatArg(FormatExpression._priority_prefix, expr.args[0])]), FormatExpression._priority_prefix
        elif isinstance(expr.function, PostfixUnaryOperator):
            return entangleLists([self.formatArg(FormatExpression._priority_postfix, expr.args[0]), [expr.function.name]]), FormatExpression._priority_postfix
        elif isinstance(expr.function, BinaryOperator):
            priority = FormatExpression._priority_binary[expr.function.name]
            args = [self.formatArg(priority, arg) for arg in expr.args]
            return entangleLists([args[0], [' ' + expr.function.name + ' '], args[1]]), priority
        elif isinstance(expr.function, BracketOperator):
            arg = self.formatArg(FormatExpression._priority_postfix, expr.args[0])
            subscript = self.formatArg(FormatExpression._priority_none, expr.args[1])
            return entangleLists([arg, ['[ '], subscript, [' ]']]), FormatExpression._priority_postfix
        else:
            raise Exception('Unknown operator: ' + repr(expr.function))

    def this(self, expr):
        return ['this'], FormatExpression._priority_terminal

    def variable(self, expr):
        return [expr.name], FormatExpression._priority_terminal

    def unformatted(self, expr):
        return [expr.value], FormatExpression._priority_none


def formatExpression(expr):
    lines, _ = FormatExpression()(expr)
    return lines
