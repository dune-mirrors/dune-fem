# %%
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
import numpy, ufl
from dune.alugrid import aluConformGrid as view
from dune.grid import cartesianDomain, Marker
from dune.fem.function import gridFunction
from dune.fem.space import lagrange, dgonb
from dune.fem.view import geometryGridView, adaptiveLeafGridView
import dune.fem
from dune.ufl import Constant
import dune.common.pickle
from transform import exact

def test1(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[10,10]) )
    grid = adaptiveLeafGridView( grid )

    """
    # mark further down fails on geometryGV due to wrong entity type from hierarchicalGrid
    x = ufl.SpatialCoordinate(ufl.triangle)
    expr =  [ (x[0]+x[1])/ufl.sqrt(2), (-x[0]+x[1])*ufl.sqrt(2) ]
    coord = lagrange(grid,dimRange=2).interpolate(expr,name="coord")
    grid = geometryGridView( coord )
    """

    @gridFunction(grid, name="gf", order=3)
    def gf(x): return numpy.sqrt(1-(x[0]**2+x[1]**2)) if x[0]**2+x[1]**2<1 else 0

    space = lagrange(grid, order=2)
    df = space.interpolate(gf,name="test")

    dune.fem.globalRefine(4,grid.hierarchicalGrid)
    # grid.hierarchicalGrid.globalRefine(4)
    for i in range(5):
        grid.hierarchicalGrid.mark(lambda e:
             Marker.refine if df.localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 1
             else Marker.coarsen)
        dune.fem.adapt(grid.hierarchicalGrid)
        df.interpolate( gf )
    print("size of adapted grid:", grid.size(0))
    # df.plot()

    # there is an issue with GeometryGV and adaptivity - perhaps
    # one needs to change the order, i.e., can only use GeometryGV<Adaptiv>
    # not the other way around?
    """
    x = ufl.SpatialCoordinate(ufl.triangle)
    expr =  [ (x[0]+x[1])/ufl.sqrt(2), (-x[0]+x[1])*ufl.sqrt(2) ]
    coord = lagrange(grid,dimRange=2).interpolate(expr,name="coord")
    grid = geometryGridView( coord )
    """

    t = Constant(0,"time")
    x   = ufl.SpatialCoordinate(ufl.triangle)
    lag = lagrange(grid, order=3)
    # lag = dgonb(grid, order=3)
    dg  = dgonb(grid, order=4, dimRange=2)
    dg5  = dgonb(grid, order=4, dimRange=4)
    df2 = lag.interpolate(exact(grid),name="exact_h")
    df5 = dg5.interpolate( [ufl.dot(x,x), -x[1],x[0], ufl.sin(x[0]*x[1])], name="euler")
    df3 = dg.interpolate(ufl.tanh(2*(t*x[0]-x[1]))*x,name="tanh")

    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump([1,2,df2,3,df3,df5],f) # adding some numbers just for testing

    series = dune.common.pickle.SeriesPickler(fileName, [1,2,df2,3,df3,df5])

    tsps = [0,0.1,0.4,0.8,2,4]
    for i,tsp in enumerate(tsps):
        t.value = tsp
        df2.interpolate(exact(grid,tsp))
        df3.interpolate(ufl.tanh(2*(t*x[0]-x[1]))*x)
        series.dump({"time":tsp})

test1("dump")
