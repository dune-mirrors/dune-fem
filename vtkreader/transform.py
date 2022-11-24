from dune.fem.function import uflFunction, gridFunction
from ufl import *
import numpy

def gradx(gv,dfs):
    g = uflFunction(gv,name="grad",order=1,ufl=grad(dfs[0])[0])
    return [g]
def grady(gv,dfs):
    g = uflFunction(gv,name="grad",order=1,ufl=grad(dfs[0])[1])
    return [g]

def exact(grid):
    @gridFunction(grid, name="exact", order=3)
    def _exact(x):
        return numpy.sin(x[0]*x[1]*numpy.pi)
    return _exact
def error(gv,dfs):
    err = uflFunction(gv,name="error",order=1,ufl=dfs[0]-exact(gv))
    return [err,*dfs]
register = [error,gradx,grady]
