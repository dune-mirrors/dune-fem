import sys
import ufl
from dune.alugrid import aluConformGrid as view
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
    print( grid.hierarchicalGrid.refineStepsForHalf )

    #####################################################
    # we can remove either of these lines to make it work
    #####################################################
    # grid = adaptiveLeafGridView( grid )              #1
    globalRefine(1,grid.hierarchicalGrid)              #2
    grid.plot()
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
    globalRefine(1,grid.hierarchicalGrid)              #2
    grid.plot()

test1("testAGV.dbf")

with open("testAGV.dbf","rb") as f:
    lag1,lag3,dg = dune.common.pickle.load(f)
grid = lag1.gridView
print( grid.hierarchicalGrid.refineStepsForHalf )
grid.plot()
elem = grid.elements.__next__()
indset = grid.indexSet
for elem in grid.elements:
    print("Elem:",elem.geometry.center,end="\t")
    print("IS:", [ indset.subIndices(elem,c) for c in range(3) ] )

# ISSUE: refinening here does not work
# RuntimeError: NotImplemented [numBlocks:/home/dedner/DUNE_PYTHON/dune-fem/dune/fem/space/mapper/codimensionmapper.hh:202]: Method numBlocks() called on non-adaptive block mapper
# globalRefine(1,grid.hierarchicalGrid)
# dune.fem.globalRefine(1,grid.hierarchicalGrid)
# grid.plot()

"""
# level = 0 works fine as expected
level=0                      #3
lag1.plot(level=level)
lag3.plot(level=level)       # this one fails
dg.plot(level=level)

globalRefine(1,[lag1,lag3,dg])
lag3.plot()
"""
