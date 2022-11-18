import dune.common
from dune.grid import structuredGrid
from dune.fem.space import lagrange, dgonb
import pickle,sys

# avoid issue with import path dune.generated missing
x = dune.common.FieldVector([1])
# avoid issue with typeRegistry for now

if sys.argv[1] == 'dump':
    grid = structuredGrid([0,0],[1,1],[10,10])
    spcL = lagrange(grid,order=4)
    spcD = dgonb(grid,order=2)
    l_h = spcL.interpolate(0,name="lag")
    d_h = spcD.interpolate(0,name="lag")
    with open("dump","wb") as f:
        # pickle.dump(grid,f)
        pickle.dump([l_h.__impl__,d_h.__impl__],f)
else:
    with open("dump","rb") as f:
        # gv = pickle.load(f)
        l_h,d_h = pickle.load(f)
print("====================")
print(l_h.size)
print(l_h.space.size)
print(l_h.space.gridView.size(0))
print(dir(l_h))
l_h.plot()
print("#######################")
