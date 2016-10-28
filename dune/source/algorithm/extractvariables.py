from __future__ import division, print_function, unicode_literals

from ..expression import Application, ConstantExpression, ConstructExpression, Expression, InitializerList, LambdaExpression, Variable


def extractVariables(expr):
    if isinstance(expr, (list, tuple)):
        return set.union(*[extractVariables(arg) for arg in expr])
    if isinstance(expr, Expression):
        if isinstance(expr, Application):
            if expr.args is not None:
                return set.union(*[extractVariables(arg) for arg in expr.args])
            else:
                return {}
        elif isinstance(expr, ConstantExpression):
            return {}
        elif isinstance(expr, ConstructExpression):
            if expr.args is not None:
                return set.union(*[extractVariables(arg) for arg in expr.args])
            else:
                return {}
        elif isinstance(expr, InitializerList):
            return set.union(*[extractVariables(arg) for arg in expr.args])
        elif isinstance(expr, LambdaExpression):
            raise Exception("Extracting used variables from lambda expression is not supported, yet.")
        elif isinstance(expr, Variable):
            return {expr}
    else:
        raise Exception("Extracting used variables is currently only supported for expressions.")
