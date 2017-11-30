from __future__ import division, print_function, unicode_literals

from ..expression import Application, ConditionalExpression, ConstantExpression, ConstructExpression, Expression, InitializerList, LambdaExpression, NullPtr, This, UnformattedExpression, Variable
from ..cplusplus import Declaration, ReturnStatement, SwitchStatement, Using


def extractVariablesFromExpression(expr):
    if isinstance(expr, Application):
        if expr.args is not None:
            return set.union(*[extractVariablesFromExpression(arg) for arg in expr.args])
        else:
            return set()
    elif isinstance(expr, ConditionalExpression):
        return set.union(*[extractVariablesFromExpression(e) for e in (expr.cond, expr.true, expr.false)])
    elif isinstance(expr, ConstantExpression):
        return set()
    elif isinstance(expr, ConstructExpression):
        if expr.args is not None:
            return set.union(*[extractVariablesFromExpression(arg) for arg in expr.args])
        else:
            return set()
    elif isinstance(expr, InitializerList):
        return set.union(*[extractVariablesFromExpression(arg) for arg in expr.args])
    elif isinstance(expr, LambdaExpression):
        return set(expr.capture)
    elif isinstance(expr, NullPtr):
        return set()
    elif isinstance(expr, This):
        return set()
    elif isinstance(expr, Variable):
        return {expr}
    elif isinstance(expr, UnformattedExpression):
        return set(expr.uses) if expr.uses is not None else set()
    else:
        raise Exception("Unknown expression", expr)


def extractVariablesFromExpressions(exprs):
    return set.union(*[extractVariablesFromExpression(expr) for expr in exprs])


def extractVariablesFromStatement(stmt):
    if isinstance(stmt, Expression):
        return extractVariablesFromExpression(stmt)
    elif isinstance(stmt, Declaration):
        if not isinstance(stmt.obj, Variable):
            raise Exception('Only variables can be declared for now.')
        return extractVariablesFromExpression(stmt.initializer) if stmt.initializer is not None else set()
    elif isinstance(stmt, ReturnStatement):
        return extractVariablesFromExpression(stmt.expression);
    elif isinstance(stmt, SwitchStatement):
        result = {stmt.var}
        for case, code in stmt.branches.items():
            result = result | extractVariablesFromStatements(code.content)
        result = result | extractVariablesFromStatements(stmt.default.content)
        return result
    elif isinstance(stmt, Using):
        return set()
    else:
        raise Exception("Unknown statement:", stmt)


def extractVariablesFromStatements(stmts):
    declared = set()
    used = set()
    for stmt in stmts:
        used = used | extractVariablesFromStatement(stmt)
        if isinstance(stmt, Declaration):
            if stmt.obj in declared:
                raise Exception("Multiply declared variable: " + stmt.obj.name)
            if stmt.obj in used:
                raise Exception("Variable declared after use: " + stmt.obj.name)
            declared.add(stmt.obj)
    return used
