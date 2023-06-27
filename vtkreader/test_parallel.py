import dune.common.pickle

from mpi4py import MPI
rank = MPI.COMM_WORLD.Get_rank()

# from dune.alugrid import aluCubeGrid as view
from dune.alugrid import aluConformGrid as view
# from dune.alugrid import aluSimplexGrid as view
from dune.grid import cartesianDomain

def dump(fileName):
    fname = f"{fileName}.{rank}.pdbf"
    grid = view( cartesianDomain([-2,-2],[2,2],[10,10]) )
    # grid.hierarchicalGrid.globalRefine(1)
    grid.hierarchicalGrid.loadBalance()
    with open(fname,"wb") as f:
        dune.common.pickle.dump([grid],f)
    return fname

print("=====================\n OUTPUT",flush=True)
fname = dump("testpara")

print("=====================\n INPUT",flush=True)
with open(fname,"rb") as f:
    dump = dune.common.pickle.load(f)
print(dump)
dump[0].plot()
