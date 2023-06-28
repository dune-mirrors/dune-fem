import sys
import ufl
from dune.alugrid import aluSimplexGrid as view
from dune.grid import cartesianDomain
from dune.fem import globalRefine
from dune.fem.space import lagrange, dgonb
from dune.fem.view import adaptiveLeafGridView
import dune.common.pickle

# issue with higher order lagrange on a globally refined adaptive LGV
# Works fine with the same space with no refinement is done (#2) and
# works fine with the same space without adaptive LGV (#1)
# Also works fine if no subsampling is used (#3)

def test1(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[2,2]) )
    print(grid.size(0))
    print(grid._module)

    #####################################################
    # we can remove either of these lines to make it work
    #####################################################
    grid = adaptiveLeafGridView( grid )                #1
    print(grid._module)
    print(grid.size(0))
    globalRefine(1,grid.hierarchicalGrid)              #2
    print(grid.size(0))
    #####################################################

    x    = ufl.SpatialCoordinate(ufl.triangle)
    lag1 = lagrange(grid, order=1)
    lag3 = lagrange(grid, order=2)                      # this one is a problem
    dg   = dgonb(grid, order=4)
    df_lag1 = lag1.interpolate(x[0]*x[1],name="lag1")
    df_lag3 = lag3.interpolate(x[0]*x[1],name="lag3")   # this one is a problem
    df_dg   =   dg.interpolate(x[0]*x[1],name="dg")

    with open(fileName,"wb") as f:
        dune.common.pickle.dump([df_lag1,df_lag3,df_dg],f)

    # globalRefine(1,[df_lag1,df_lag3,df_dg])
    # df_lag3.plot(level=1)
    indset = grid.indexSet
    for elem in grid.elements:
        print("Elem:",elem.geometry.center,end="\t")
        print("IS:", [ indset.subIndices(elem,c) for c in range(3) ] )

test1("testAGV.dbf")

with open("testAGV.dbf","rb") as f:
    lag1,lag3,dg = dune.common.pickle.load(f)
grid = lag1.gridView
elem = grid.elements.__next__()
indset = grid.indexSet
for elem in grid.elements:
    print("Elem:",elem.geometry.center,end="\t")
    print("IS:", [ indset.subIndices(elem,c) for c in range(3) ] )

"""
# level = 0 works fine as expected
level=0                      #3
lag1.plot(level=level)
lag3.plot(level=level)       # this one fails
dg.plot(level=level)

globalRefine(1,[lag1,lag3,dg])
lag3.plot()
"""
