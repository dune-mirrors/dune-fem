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

    grid.hierarchicalGrid.globalRefine(4)
    for i in range(5):
        grid.hierarchicalGrid.mark(lambda e:
             Marker.refine if df.localFunction(e).jacobian([1./3.,1./3.]).infinity_norm > 1
             else Marker.coarsen)
        dune.fem.adapt(grid.hierarchicalGrid)
        df.interpolate( gf )
    print("size of adapted grid:", grid.size(0))
    # df.plot()

    x = ufl.SpatialCoordinate(ufl.triangle)
    expr =  [ (x[0]+x[1])/ufl.sqrt(2), (-x[0]+x[1])*ufl.sqrt(2) ]
    coord = lagrange(grid,dimRange=2).interpolate(expr,name="coord")
    grid = geometryGridView( coord )

    dg = dgonb(grid, order=4)
    df2 = dg.interpolate(exact(grid),name="test2")
    # df2.plot()

    with open(fileName,"wb") as f:
        dune.common.pickle.dump([1,2,df2,3],f) # adding some numbers just for testing
    with open("grid"+fileName,"wb") as f:
        dune.common.pickle.dump([grid],f) # adding some numbers just for testing
    df2.gridView.writeVTK("dump", pointdata=[df2])
    df2.gridView.writeVTK("dump2", pointdata=[df2], subsampling=2)

test1("dump.dbf")
