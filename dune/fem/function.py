import math

class MathExpression(object):
    def __init__(self, expr):
        self.expr = expr
    def __call__(self,x):
        dimD = len(x)
        x0 = x[0]
        x1 = 0 if dimD < 2 else x[1]
        x2 = 0 if dimD < 3 else x[2]
        y = []
        for e in self.expr:
            y += [eval(e)]
        return y
