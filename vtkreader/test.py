import sys,math

# from dunereader import setDuneModulePaths
# setDuneModulePaths()

# gv = "geometry"
# gv = "adaptive"
gv = "original"
def gridView(grid):
    if gv == "geometry":
        expr =  [ (x[0]+x[1])/ufl.sqrt(2), (-x[0]+x[1])*ufl.sqrt(2) ]
        # coord = uflFunction(grid,name="coord",order=1,ufl=expr)
        coord = lagrange(grid,dimRange=2).interpolate(expr,name="coord")
        return geometryGridView(coord)
    elif gv == "adaptive":
        return adaptiveLeafGridView(grid)
    else:
        return grid

if sys.argv[1] == 'dump':
    import ufl
    x = ufl.SpatialCoordinate(ufl.triangle)

    import dune.common.pickle
    from dune.grid import structuredGrid
    from dune.fem.space import lagrange, dgonb
    from dune.fem.function import uflFunction
    from dune.fem.view import geometryGridView, adaptiveLeafGridView
    grid = gridView( structuredGrid([0,0],[1,1],[10,10]) )
    spcL = lagrange(grid,order=4)
    spcD = dgonb(grid,order=2)
    l_h = spcL.interpolate(x[0]*x[1],name="lag")
    d_h = spcD.interpolate(ufl.tanh(5*x[0]*x[1]),name="lag")
    with open("dump"+gv,"wb") as f:
        dune.common.pickle.dump([1,l_h,2,d_h,3],f)
else:
    import dune.common.pickle
    with open("dump"+gv,"rb") as f:
        _,l_h,_,d_h,_ = dune.common.pickle.load(f)
print("====================")
# print(dir(l_h.space.gridView))
print(l_h.size,l_h.__impl__)
print(l_h.space.size)
print(l_h.space.gridView.size(0))
# l_h.plot()
print("#######################")
