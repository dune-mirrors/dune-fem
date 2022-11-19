import dune.common
from dune.grid import structuredGrid
from dune.fem.space import lagrange, dgonb
import dune.fem
import sys

# avoid issue with import path dune.generated missing
# x = dune.common.FieldVector([1])
# avoid issue with typeRegistry for now

if sys.argv[1] == 'dump':
    grid = structuredGrid([0,0],[1,1],[10,10])
    spcL = lagrange(grid,order=4)
    spcD = dgonb(grid,order=2)
    l_h = spcL.interpolate(0,name="lag")
    d_h = spcD.interpolate(0,name="lag")
    with open("dump","wb") as f:
        dune.fem.dump([1,l_h,2,d_h,3],f)
else:
    with open("dump","rb") as f:
        _,l_h,_,d_h,_ = dune.fem.load(f)
print("====================")
print(l_h.size,l_h.__impl__)
print(l_h.space.size)
print(l_h.space.gridView.size(0))
l_h.plot()
print("#######################")
