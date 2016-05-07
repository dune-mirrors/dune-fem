import math
try:
    from numbae import jit
except ImportError:
    pass

try:
    import sympy
    class SympyExpression(object):
        def __init__(self, expr):
            self.dimR = len(expr)
            try:
                self.expr_func = jit()(SympyExpression.evaluate_)
                print( "using numba.jit" )
            except NameError:
                self.expr_func = SympyExpression.evaluate_
            self.exp_sympy = [None]*self.dimR
            x0, x1, x2 = sympy.symbols('x0 x1 x2')
            for r in range(0,self.dimR):
                self.exp_sympy[r] = eval(expr[r])
        def evaluate(self,x,res):
            self.expr_func(self,x,res)
        def evaluate_(self,x,res):
            x0, x1, x2 = sympy.symbols('x0 x1 x2')
            dimD = len(x)
            for r in range(0,self.dimR):
                resExp = self.exp_sympy[r]
                for i in range(0,dimD):
                    resExp = resExp.subs("x"+str(i),x[i])
                res[r] = float( resExp.evalf() )
except ImportError:
    pass

class MathExpression(object):
    def __init__(self, expr):
        self.expr = expr
        self.dimR = len(expr)
    def evaluate(self,x,res):
        dimD = len(x)
        x0 = x[0]
        if (dimD>1):
            x1 = x[1]
        else:
            x1 = 0
        if (dimD>2):
            x2 = x[2]
        else:
            x2 = 0
        for r in range(0,self.dimR):
            res[r] = eval(self.expr[r])
class FuncExpression(object):
    def __init__(self, dimR, expr):
        self.dimR = dimR
        try:
            self.expr = jit()(expr)
            print( "using numba.jit" )
        except NameError:
            self.expr = expr
    def evaluate(self,x,res):
        self.expr(x,res)
class LocalFuncExpression(object):
    def __init__(self, dimR, expr):
        self.dimR = dimR
        try:
            self.expr = jit()(expr)
            print( "using numba.jit" )
        except NameError:
            self.expr = expr
    def evaluate(self,entity,x,res):
        self.expr(entity,x,res)
