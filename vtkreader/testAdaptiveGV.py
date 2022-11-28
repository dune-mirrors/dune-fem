import ufl
from dune.alugrid import aluConformGrid as view
from dune.grid import cartesianDomain
from dune.fem import globalRefine
from dune.fem.space import lagrange, dgonb
from dune.fem.view import adaptiveLeafGridView
import dune.common.pickle

def test1(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[10,10]) )

    #####################################################
    # we can remove either of these lines to make it work
    #####################################################
    grid = adaptiveLeafGridView( grid )
    print(grid.size(0))
    globalRefine(1,grid.hierarchicalGrid)
    print(grid.size(0))
    # grid.hierarchicalGrid.globalRefine(1)
    #####################################################

    x    = ufl.SpatialCoordinate(ufl.triangle)
    lag1 = lagrange(grid, order=1)
    lag3 = lagrange(grid, order=3)
    dg   = dgonb(grid, order=4)
    df_lag1 = lag1.interpolate(x[0]*x[1],name="lag1")
    df_lag3 = lag3.interpolate(x[0]*x[1],name="lag3")
    df_dg   =   dg.interpolate(x[0]*x[1],name="dg")

    with open(fileName,"wb") as f:
        dune.common.pickle.dump([df_lag1,df_lag3,df_dg],f)

test1("testAGV.dbf")


with open("testAGV.dbf","rb") as f:
    lag1,lag3,dg = dune.common.pickle.load(f)

# level = 0 works fine as expected
level=1
lag1.plot(level=level)
lag3.plot(level=level) # this one fails
dg.plot(level=level)
