import io, pickle, sys
from mpi4py import MPI
import sys
try:
    import petsc4py
    petsc4py.init(sys.argv)
    from petsc4py import PETSc
    print("""***************************************************
*** parallelTest.py: petsc4py enabled!
***************************************************
          """)
except:
    petsc4py = None
    print("""***************************************************
***   WARNING:
***
***   parallelTest.py: petsc4py was disabled due to ImportError!
***************************************************
          """)

assert petsc4py is not None, "test whether petsc4py is really imported!"

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

from dune.commands import makegenerated
from dune.generator import algorithm, requiredModules
from dune.grid import yaspGrid, cartesianDomain, Partitions
from dune.alugrid import aluConformGrid as confGrid
from dune.alugrid import aluSimplexGrid as simpGrid
from dune.alugrid import aluCubeGrid as cubeGrid
from dune.fem import assemble
from dune.fem.plotting import plotPointData as plot
from dune.fem.space import combined
from dune.fem.scheme import galerkin as solutionScheme
from dune.fem.view import adaptiveLeafGridView as adaptive
from dune.fem.view import filteredGridView as filterGV
from dune.fem import integrate
from dune.ufl import Constant, DirichletBC
from ufl import ( TestFunction, TrialFunction, SpatialCoordinate, triangle, FacetNormal,
                  dx, ds, grad, div, grad, dot, inner, sqrt, exp, conditional,
                  as_vector, avg, jump, dS, CellVolume, FacetArea, atan )

def compute(comm,useAdapt,gridType,
            name, dirichlet, dg, storage,
            useCombined=False,
            internal=True, external=(petsc4py is not None),
            verbose=True,
            showPlot=False, ignore=False):

    global Print
    if not verbose:
        Print = emptyPrint

    iterations = []
    failure = False

    testName = ("adapt" if useAdapt else "" ) + f"{gridType}-{name}"
    # outFile = io.StringIO()
    outFile = sys.stdout

    if gridType == "conf" and dg: return [],[]               # this can not work due to missing ghosts

    Print(comm,"--------------------------", file=outFile)
    Print(comm,f"Running: {testName}", file=outFile)
    if ignore:
        Print(comm,"Ignored...", file=outFile)
        Print(comm,"--------------------------", file=outFile)
        return

    domain = cartesianDomain( [-1, -1], [1, 1], [21, 21], overlap=1)
    gridName = gridType + "Grid"
    # create grid view, i.e. yaspGrid, cubeGrid, confGrid, simpGrid
    gridView = eval(gridName)(domain)

    if useAdapt:
        gridView = adaptive( gridView ) # , partition=Partitions.interiorBorder )

    if dg:
        from dune.fem.space import dgonb as Space # lagrange as Space
    else:
        from dune.fem.space import lagrange as Space

    if not useCombined:
        dimR = 3
        space = Space(gridView, order=3, dimRange=dimR, storage=storage)
    else:
        dimR = 3
        spc1 = Space(gridView, dimRange=2, order=3, storage=storage)
        spc2 = Space(gridView, order=2, storage=storage)
        space = combined(spc1,spc2)

    # %%
    x    = SpatialCoordinate(space)
    u    = TrialFunction(space)
    v    = TestFunction(space)
    n    = FacetNormal(space)
    he   = avg( CellVolume(space) ) / FacetArea(space)
    hbnd = CellVolume(space) / FacetArea(space)

    # diffusion factor
    eps = Constant(0.01,"eps")
    # dg and weak boundary penalty
    beta = Constant(10 * space.order**2 if space.order > 0 else 1, name="beta")

    exact = atan(10*x[1]*(1-x[0])**2/4)
    forcing = -div(eps*grad(exact)) + exact
    g = exact

    aInternal     = dot(eps*grad(u[0]), grad(v[0])) * dx
    diffSkeleton  = eps*beta/he*jump(u[0])*jump(v[0])*dS -\
                    eps*dot(avg(grad(u[0])),n('+'))*jump(v[0])*dS -\
                    eps*jump(u[0])*dot(avg(grad(v[0])),n('+'))*dS
    diffBnd       = eps*beta/hbnd*(u[0]-g)*v[0]*ds -\
                    eps*dot(grad(u[0]),n)*v[0]*ds -\
                    eps*(u[0]-g)*dot(grad(v[0]),n)*ds

    form = aInternal

    for r in range(0,dimR):
        form += inner(u[r],v[r]) * dx

    if dg:
        assert not dirichlet
        form += diffSkeleton

    if dirichlet:
        dbcs = [ DirichletBC(space, [g]+(dimR-1)*[0], b) for b in [1,2,3,4] ]
    else:
        form += diffBnd
        dbcs  = []

    eqn = (form == forcing*v[0]*dx)
    scheme = solutionScheme([eqn, *dbcs], solver="cg",
                parameters={"nonlinear.verbose":False,
                            "nonlinear.maxiterations":10,
                            "linear.preconditioning.method":"jacobi",
                            "linear.verbose":False,
                            "linear.tolerance":1e-10,
                            "linear.gmres.restart":100,
                            "linear.maxiterations":5000})

    uh1 = space.interpolate(dimR*[14], name="solution")

    # check storage presence
    assert hasattr( uh1, "as_"+storage ), f"as_{storage} should be present if storage {storage} was chosen!"

    code = "template <class DF> void finalize(DF &df) { df.communicate(); }"
    communicate = algorithm.load('finalize', io.StringIO(code), uh1)

    rng = np.random.default_rng()
    if storage == "numpy":
        dofSize = uh1.as_numpy[:].shape[0]
    elif external and storage == "petsc":
        try:
            dofSize = uh1.as_petsc.getArray().shape[0]
        except AttributeError: # ci gives numpy.ndarray object has no attribute getArray
            dofSize = uh1.as_petsc.shape # not sure why this is right?
    else:
        dofSize = -1
    Print(comm, "SIZES:", space.size, uh1.size, dofSize, all=True)
    if internal:
        Print(comm,"internal dune-fem solver", file=outFile)
        if storage == "numpy":
            uh1.as_numpy[:] = rng.standard_normal(dofSize)
            communicate(uh1)
        """
        # this does not work as expected
        if storage == "petsc":
            uh1.as_petsc.getArray()[:] = rng.standard_normal(dofSize)
            communicate(uh1)
        """
        try:
            info = scheme.solve(target=uh1)
            Print(comm, "solver internal:",info["converged"], info["linear_iterations"], info["iterations"], file=outFile)
            if showPlot: uh1.plot(clim=[-1.5,1.5], level=3, partition=Partitions.all)
            error1 = np.sqrt( integrate((uh1[0]-exact)**2, order=6 ))
            Print(comm, "error:",error1, all=True)
            info1 = scheme.solve(target=uh1)
            if info1["linear_iterations"] > 1:
                print("Error: second internal solver required",info1["linear_iterations"],"iterations")
                info["converged"] = False
            iterations += [(testName,info["linear_iterations"])]
        except:
            info = {'converged':False}
    else:
        #assert False
        #sys.exit(0)
        Print(comm,"internal not tested", file=outFile)
        Print(comm,"===================", file=outFile)
        error1 = 0
        info = {'converged':True, "iterations":1}
    if external and storage == "petsc":
        # b = space.function(name="rhs")
        # A = scheme.linear()
        # scheme.jacobian(space.zero, A,b)
        A,b = assemble([eqn,*dbcs])
        Amat = A.as_petsc

        ksp = PETSc.KSP().create()
        ksp.setType(PETSc.KSP.Type.CG)
        pc = ksp.getPC()
        pc.setType("jacobi")  # default seems to be ilu and that fails with lagrange+Dirichlet but 'none' fails with weak Dirichlet
        ksp.setOperators(Amat, Amat)
        ksp.setTolerances(rtol=1e-7,atol=1e-10)
        ksp.setFromOptions()

        def petscSolver():
            ksp.setConvergenceHistory()
            ksp.setInitialGuessNonzero(True)
            if dirichlet: # I don't understand why this fails with weak Dirichlet
                scheme.setConstraints(uh2)
            # not sure why these tests work without the communications
            continuous = False
            if continuous: # need to communicate on continuous spaces
                code = "template <class DF> void prepare(DF &df) { df.dofVector().clearGhost(); }"
                clearGhost = algorithm.load('prepare', io.StringIO(code), uh1)
                clearGhost(uh2)
            ksp.solve(b.as_petsc,uh2.as_petsc)
            if continuous: # need to communicate on continuous spaces
                communicate(uh2)
            return len(ksp.getConvergenceHistory())


        uh2 = space.interpolate(dimR*[14], name="solution")
        # this does not work as expected
        #uh2.as_petsc.getArray()[:] = rng.standard_normal(dofSize)
        #communicate(uh2)
        petscIter = petscSolver()
        if showPlot: uh2.plot(allowNaN=True, clim=[-1.5,1.5])
        # pyplot.semilogy(ksp.getConvergenceHistory())
        Print(comm,"solver petsc4py:", petscIter, file=outFile)
        petscIter1 = petscSolver()
        if petscIter1 > 1:
            print("Second solver required",petscIter1,"iterations")
            uh2.clear() # this will mark as failed
        iterations += [(testName+"-petsc4py",info["linear_iterations"])]
    else:
        petscIter = None
        Print(comm,"external not tested", file=outFile)
        # assert False

    if external and internal and storage == "petsc":
        if showPlot: plot(uh1-uh2, grid=gridView)
        error = uh1.as_petsc.array - uh2.as_petsc.array
        difference = comm.reduce(np.dot(error,error), op=MPI.SUM)
        Print(comm,"Difference", difference, file=outFile)
    else:
        difference = 0

    Print(comm,"--------------------------", file=outFile)

    # we would like to compare this between different number of procs
    fail = 0
    fail +=   1 if not info["converged"] else 0
    fail +=  10 if error1 > 1e-4 else 0
    fail += 100 if info["iterations"]>1 else 0
    if petscIter is not None:
        fail +=  1000 if petscIter > 600 else 0
    if difference is not None:
        fail += 10000 if difference > 1e-5 else 0
    if fail > 0:
        if not outFile == sys.stdout:
            Print(comm, outFile.getvalue())
        Print(comm,"--------------------------")
        Print(comm,f"{testName}: Failed, reason {fail}")
        Print(comm,"--------------------------")
        failures = [(testName,fail)]
    else:
        Print(comm,"--------------------------")
        Print(comm,f"{testName}: Success")
        Print(comm,"--------------------------")
        failures = []
    return failures, iterations

############################################################################
def experiments(ignore):
    return [
      # can add things like:
      # "external":True, "showPlot":False, "ignore":False},
      { "name":"dg-petsc",
        "dg":True, "dirichlet":False, "storage":"petsc",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lag-petsc",
        "dg":False, "dirichlet":False, "storage":"petsc",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lagDir-petsc",
        "dg":False, "dirichlet":True, "storage":"petsc",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"dg-numpy",
        "dg":True, "dirichlet":False, "storage":"numpy",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lag-numpy",
        "dg":False, "dirichlet":False, "storage":"numpy",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lagDir-numpy",
        "dg":False, "dirichlet":True, "storage":"numpy",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"dg-istl",
        "dg":True, "dirichlet":False, "storage":"istl",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lag-istl",
        "dg":False, "dirichlet":False, "storage":"istl",
        "showPlot":False, "ignore":ignore,
      },
      { "name":"lagDir-istl",
        "dg":False, "dirichlet":True, "storage":"istl",
        "showPlot":False, "ignore":ignore,
      },
    ]


def runTest(count=None):
    if count is None:
        count = datetime.date.today().toordinal()
    grids = [ ["yasp","cube"], ["conf","simp"]]
    storage = [ "numpy", "petsc", "istl" ]
    crossProduct = [ [g,s] for g in grids for s in storage ]

    # we have 12 different sets of experiments currently so everything is
    # checked in about 2 weeks
    total = len(crossProduct)*2  # the 2 is for with or without 'adaptive leaf grid view"
    todaysExperiment = count % total
    useAdapt  = todaysExperiment % 2
    gridGroup = crossProduct[todaysExperiment // 2][0]
    storage   = crossProduct[todaysExperiment // 2][1]

    comm = MPI.COMM_WORLD
    ex = experiments(ignore=False)
    ex = [e for e in ex if e["storage"] == storage]

    # in serial first compile all required modules using multiple threads without computing any results
    if comm.Get_size() == 1:
        with ThreadPoolExecutor(max_workers=4) as executor:
            for e in ex:
                for g in gridGroup:
                    executor.submit(compute, comm, gridType=g, useAdapt=useAdapt, **e,
                                    internal=False, external=False,verbose=False)

    # now the actual computation
    iterations = []
    failures = []
    for e in ex:
        for g in gridGroup:
            f,i = compute(comm, gridType=g,useAdapt=useAdapt, **e)
            failures += f
            iterations += i

    Print(comm, failures)
    Print(comm, iterations)

    mpiSize = comm.Get_size()
    with open(f"parallelTest{mpiSize}.out","bw") as f:
        pickle.dump([failures,iterations],f)

    return len(failures)==0


def main():
    if len(sys.argv) > 1:
        count = int(sys.argv[1])
        return runTest(count)
    iterations = []
    failures = []
    grids = ["yasp","conf","simp","cube"]
    useAdapt = True
    if useAdapt:
        makegenerated(fileName="modulesAdapt.txt", threads=12)
    else:
        makegenerated(fileName="moduleslist.txt", threads=12)

    comm = MPI.COMM_WORLD
    ignore = False
    ex = experiments(ignore)
    for e in ex:
        Print(comm,"******************************************************************")
        for g in grids:
            f,i = compute(comm,gridType=g,useAdapt=useAdapt, **e)
            failures += f
            iterations += i
        Print(comm,"******************************************************************")

    print(iterations)
    print(failures)

    Print(comm,"*** main finished! ***")
    requiredModules("tmp.txt")
    if comm.Get_rank() == 0:
        if useAdapt:
            requiredModules("modulesAdapt.txt") # outcomment once generated in case you change the 'ignore' values above
            pass
        else:
            requiredModules("moduleslist.txt")
            pass

if __name__ == "__main__":
    print(sys.argv)
    if main():
        sys.exit(0)
    else:
        sys.exit(1)
