from __future__ import division, print_function, unicode_literals

from ..builtin import BuiltInFunction
from ..expression import Application, ConditionalExpression, ConstantExpression, ConstructExpression, Expression, InitializerList, LambdaExpression, UnformattedExpression, Variable
from ..cplusplus import Declaration, ReturnStatement, Using


def extractIncludesFromExpression(expr):
    if isinstance(expr, Application):
        if isinstance(expr.function, BuiltInFunction):
            result = {expr.function.header}
        else:
            result = set()
        if expr.args is not None:
            result.update(set.union(*[extractIncludesFromExpression(arg) for arg in expr.args]))
        return result
    elif isinstance(expr, ConstantExpression):
        return set()
    elif isinstance(expr, ConstructExpression):
        if expr.args is not None:
            return set.union(*[extractIncludesFromExpression(arg) for arg in expr.args])
        else:
            return set()
    elif isinstance(expr, InitializerList):
        return set.union(*[extractIncludesFromExpression(arg) for arg in expr.args])
    elif isinstance(expr, LambdaExpression):
        return set()
    elif isinstance(expr, Variable):
        return set()
    elif isinstance(expr, UnformattedExpression):
        return set()
    elif isinstance(expr, ConditionalExpression):
        return set.union(*[extractIncludesFromExpression(arg) for arg in (expr.cond, expr.true, expr.false)])
    else:
        raise Exception("Unknown expression", expr)


def extractIncludesFromExpressions(exprs):
    return set.union(*[extractIncludesFromExpression(expr) for expr in exprs]) if exprs else set()


def extractIncludesFromStatement(stmt):
    if isinstance(stmt, Expression):
        return extractIncludesFromExpression(stmt)
    elif isinstance(stmt, Declaration):
        if not isinstance(stmt.obj, Variable):
            raise Exception('Only variables can be declared for now.')
        return extractIncludesFromExpression(stmt.initializer) if stmt.initializer is not None else set()
    elif isinstance(stmt, ReturnStatement):
        return extractIncludesFromExpression(stmt.expression);
    elif isinstance(stmt, Using):
        if isinstance(stmt.obj, BuiltInFunction):
            return {stmt.obj.header}
        else:
            return set()
    else:
        raise Exception("Unknown statement:", stmt)


def extractIncludesFromStatements(stmts):
    return set.union(*[extractIncludesFromStatement(stmt) for stmt in stmts]) if stmts else set()
