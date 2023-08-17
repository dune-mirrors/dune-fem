# %%
import matplotlib
matplotlib.rc( 'image', cmap='jet' )
import numpy, ufl
from dune.alugrid import aluConformGrid as view
# from dune.alugrid import aluCubeGrid as view
from dune.grid import cartesianDomain, Marker
from dune.fem.function import gridFunction, uflFunction
from dune.fem.space import lagrange, dgonb
from dune.fem.view import geometryGridView, adaptiveLeafGridView
import dune.fem
from dune.ufl import Constant
import dune.common.pickle
from transform import exact

def adaptTest(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[4,4]) )
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

    dune.fem.globalRefine(2,grid.hierarchicalGrid)
    # grid.hierarchicalGrid.globalRefine(4)
    for i in range(5):
        grid.hierarchicalGrid.mark(lambda e:
                  Marker.refine if df.localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 1
             else Marker.coarsen)
        dune.fem.adapt(grid.hierarchicalGrid)
        df.interpolate( gf )

    t = Constant(0,"time")
    x   = ufl.SpatialCoordinate(ufl.triangle)
    lag = lagrange(grid, order=3)
    dg  = dgonb(grid, order=4, dimRange=2)
    dg5  = dgonb(grid, order=4, dimRange=4)
    df2 = lag.interpolate(exact(grid),name="exact_h")
    df5 = dg5.interpolate( [ufl.dot(x,x)/4, -x[1],x[0], ufl.sin(2*ufl.pi*x[0]*x[1])], name="euler")
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

def geoTest(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[10,10]) )
    x   = ufl.SpatialCoordinate(ufl.triangle)
    coordSpace = lagrange(grid, dimRange=2, order=4)
    coord = ufl.as_vector([ x[0]*(1+0.1*ufl.cos(ufl.pi*x[1])), x[1]*(1+0.2*ufl.sin(ufl.pi*x[0])) ])
    coord = coordSpace.interpolate(coord,name="coord")
    geoGrid = geometryGridView(coord)
    geoLag = lagrange(geoGrid, order=4)
    geoDf = geoLag.interpolate(ufl.sin(ufl.pi*ufl.dot(x,x)), name="df")
    geoDf.plot(level=4)
    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump(["hallo",geoDf],f)

def surfTest(fileName):
    grid = view( cartesianDomain([-2,-2],[2,2],[4,4]) )
    x = ufl.SpatialCoordinate(ufl.triangle)
    coord = ufl.as_vector([ x[0],x[1], ufl.sin(ufl.pi*x[0]*x[1]) ])
    coord_h = lagrange(grid, dimRange=3, order=3).interpolate(coord,name="coord")
    oord_h = uflFunction(grid, name="coord", order=5, ufl=coord)
    geoGrid = geometryGridView(coord_h)
    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump([geoGrid],f)
    geoLag = lagrange(geoGrid, order=3)
    x = ufl.SpatialCoordinate(geoLag)
    geoDf = geoLag.interpolate(ufl.sin(ufl.pi*(x[0]+x[1]+x[2])), name="df")
    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump(["hallo",geoDf,coord],f)

def test3D(fileName):
    grid = view( cartesianDomain([-2,-2,-2],[2,2,2],[4,4,4]) )
    space = lagrange(grid, order=4)
    x = ufl.SpatialCoordinate(space)
    df1 = space.interpolate(ufl.sin(2*ufl.pi*x[0]*x[1]*x[2]), name="df")
    df2 = space.interpolate(ufl.dot(x,x), name="df")
    with open(fileName+".dbf","wb") as f:
        dune.common.pickle.dump([df1,1,df2],f)

############################################################################

adaptTest("adapt")
geoTest("geo")
surfTest("surface")
test3D("3d")
