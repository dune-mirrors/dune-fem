import dune.common
from dune.grid import structuredGrid
from dune.fem.space import lagrange, dgonb
import pickle,sys

# avoid issue with import path dune.generated missing
x = dune.common.FieldVector([1])
# avoid issue with typeRegistry for now
grid = structuredGrid([0,0],[1,1],[10,10])

if sys.argv[1] == 'dump':
    spcL = lagrange(grid,order=4)
    spcD = dgonb(grid,order=2)
    with open("dump","wb") as f:
        # pickle.dump(grid,f)
        pickle.dump([spcL.__impl__,spcD.__impl__],f)
    print(grid.size(0),spcL.size,spcL.gridView.size(0),spcD.size)
else:
    with open("dump","rb") as f:
        # gv = pickle.load(f)
        spcL,spcD = pickle.load(f)
    print("====================")
    # print(gv.size(0))
    print(spcL.size,spcD.size)
    print(spcL.gridView.size(0))
    print("#######################")

    """ fails with yasp - test with alu
    gv.hierarchicalGrid.globalRefine(2)

    print("====================")
    print(gv.size(0))
    print(spc.size)
    print(spc.gridView.size(0))
    """
