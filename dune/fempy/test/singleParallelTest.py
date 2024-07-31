import io, pickle, sys
from mpi4py import MPI
import sys
try:
    import petsc4py
    petsc4py.init(sys.argv)
    from petsc4py import PETSc
except:
    petsc4py = None
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from matplotlib import pyplot

#from dune.common import comm
#print = functools.partial(print, flush=True) if comm.rank == 0 else lambda *a, **kw: None

def Print(comm,*args,**kwargs):
    allRanks = kwargs.pop("all",False)
    if comm.Get_rank() == 0 or allRanks:
        # print(f"{0}/{comm.Get_size()}",*args,**kwargs,flush=True)
        print(*args,**kwargs,flush=True)
    if allRanks:
        comm.barrier()

def emptyPrint(comm,*args,**kwargs):
    pass

from parallelTest import compute, experiments, Print

def main():
    comm = MPI.COMM_WORLD
    ex = experiments(False)
    if len(sys.argv) <= 1:
        print(f"usage: {sys.argv[0]} <grid> <experiment number>")
        print(f"example: {sys.argv[0]} yasp 1 for running adaptyasp-lag-petsc")
        ex = experiments(False)
        for i,e in enumerate(ex):
            print(f"Ex {i}: {e}")
        sys.exit(1)
        return 1

    iterations = []
    failures = []
    g = sys.argv[1] # grid
    print(g)
    e = int(sys.argv[2])

    useAdapt = True
    print(ex[e])
    compute(comm, useAdapt, g, **ex[e])


    print(iterations)
    print(failures)

    Print(comm,"*** main finished! ***")

if __name__ == "__main__":
    if main():
        sys.exit(0)
    else:
        sys.exit(1)
