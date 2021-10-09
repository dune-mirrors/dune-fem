import numpy
from dune.fem.plotting import plotPointData as plot
import dune.create as create
from dune.grid import structuredGrid
from dune.istl import blockVector
import ufl

g = structuredGrid([0,0],[1,1],[2,3])

s = create.space("lagrange",g,dimRange=2,storage="istl")
f1 = s.interpolate(expr=[2,1], name="tmp")
dofs = blockVector(int(s.size/s.localBlockSize), s.localBlockSize)
# f2 = s.function("tmp", expr=[2,1], dofVector=dofs)
f2 = s.function("tmp", dofVector=dofs)
f2.interpolate([2,1])
assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f1.as_istl)])
assert all([(d1-d2).two_norm==0 for d1,d2 in zip(dofs,f2.as_istl)])
operator = create.operator("galerkin",
                 ufl.dot(ufl.TrialFunction(s),ufl.TestFunction(s)) * ufl.dx)
f1 = s.function("tmp", [2,1], blockVector(s.size//s.localBlockSize,s.localBlockSize) )
f2 = s.function("tmp", [2,1], blockVector(s.size//s.localBlockSize,s.localBlockSize) )
operator(f1,f2)

s = create.space("lagrange",g,dimRange=2,storage="numpy")
f1 = s.interpolate([2,1], name="tmp")
dofs = numpy.ndarray(s.size)
f2 = s.function("tmp", [2,1], dofs)
assert not (dofs-f1.as_numpy).any()
assert not (dofs-f2.as_numpy).any()
f1.interpolate([3,2])
dofs[:] = f1.as_numpy
assert not (dofs-f1.as_numpy).any()
assert not (dofs-f2.as_numpy).any()

operator = create.operator("galerkin",
                 ufl.dot(ufl.TrialFunction(s),ufl.TestFunction(s)) * ufl.dx)
f1 = s.function("tmp", [2,1], numpy.ndarray(s.size))
f2 = s.function("tmp", [2,1], numpy.ndarray(s.size))
operator(f1,f2)

### the following will fail since the change in the underlying dof storage
### for f1 due to the refinement step will not be mimicked for the storage of dofs
# g.hierarchicalGrid.globalRefine(1)
# assert dofs.shape == f1.as_numpy.shape

try:
    import petsc4py
    from petsc4py import PETSc
    petscs  = create.space("lagrange",g,dimRange=2,storage="petsc")
    petscf1 = petscs.interpolate([2,1], name="tmp")
    petscDofs = PETSc.Vec().createSeq(petscs.size,bsize=2)
    # print("1.",petscDofs[0],petscDofs[1])
    petscf2 = petscs.function("tmp", [2,1], petscDofs)
    # print("2.",petscDofs[0],petscDofs[1])
    # x = ufl.SpatialCoordinate(s)
    # petscf2.interpolate([x[0]*x[1],2])
    # petscf2.plot()
    petscDiff = petscs.interpolate(petscf1-petscf2, name="tmp")
    assert petscDiff.scalarProductDofs(petscDiff) == 0
    # assert all([(d1-d2).two_norm==0 for d1,d2 in zip(petscDofs,petscf1.as_petsc)])
    # assert all([(d1-d2).two_norm==0 for d1,d2 in zip(petscDofs,petscf2.as_petsc)])
    print("NOTE: there is an issue with the tests commented out above - iteration over the petsc vectors seems broken")
except ImportError:
    pass
