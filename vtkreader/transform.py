from dune.fem.function import uflFunction, gridFunction
from ufl import *
import numpy

def gradx(gv,df,dfs):
    if df.dimRange > 1:
        g = sqrt(dot(df,df))
        name = f"d/dx |{df.name}|"
    else:
        g = df
        name = f"d/dx {df.name}|"
    g = uflFunction(gv,name=name,order=1,ufl=grad(g)[0])
    return [g]
def grady(gv,df,dfs):
    if df.dimRange > 1:
        g = sqrt(dot(df,df))
    else:
        g = df
    g = uflFunction(gv,name="grad",order=1,ufl=grad(g)[1])
    return [g]

def exact(grid):
    @gridFunction(grid, name="exact", order=3)
    def _exact(x):
        return numpy.sin(x[0]*x[1]*numpy.pi)
        # return x[0]
    return _exact
def error(gv,df,dfs):
    # err = uflFunction(gv,name="error",order=1,ufl=dfs[0]-exact(gv))
    err = uflFunction(gv,name="error",order=1,ufl=(dfs[0]-exact(gv)))
    return [err,*dfs]
def velocity(gv,df,dfs):
    U = dfs[2]
    return [uflFunction(gv,name="velo",order=1,ufl=[U[1],U[2],0])]

# this does not currently work - needs some more thought
# One issue is that the transformer is only called after the grid has been
# changed
def adapt(gv,dfs):
    from dune.grid import Marker
    from dune.fem import adapt
    for i in range(4):
        dfs[0].gridView.hierarchicalGrid.mark(lambda e:
             Marker.refine if dfs[0].localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 4
                           else Marker.keep)
        adapt(dfs)
    return dfs
register = [error,gradx,grady,adapt,velocity]
