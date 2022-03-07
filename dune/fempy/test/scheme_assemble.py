#!/usr/bin/env python3

import time,sys,io
import numpy as np
import dune.grid
import dune.fem
import dune.create as create
from dune.fem.operator import linear as linearOperator
from dune.generator import algorithm

from dune.ufl import Space
from ufl import TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad

test_numpy   = True
try:
    import dune.istl
    test_istl  = True
except ImportError:
    test_istl  = False
try:
    import petsc4py
    # test_petsc = True
    test_petsc = False # issue with petsc state
except:
    test_petsc = False

testLoop = 1
# testLoop = 10
# testLoop = 1000

# grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [89, 89]), dimgrid=2)
try:
    grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [9, 9]), dimgrid=2)
except:
    grid = create.grid("Yasp", dune.grid.cartesianDomain([0, 0], [1, 1], [9, 9]), dimgrid=2)

d = 0.001
p = 1.7

uflSpace = Space((2,2),1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())

rhs = (x[0] + x[1]) * v[0]
a = (pow(d + inner(grad(u), grad(u)), (p-2)/2)*inner(grad(u), grad(v)) + grad(u[0])[0]*v[0]) * dx + 10*inner(u, v) * ds
b = rhs * dx + 10*rhs * ds
model = create.model("integrands", grid, a==b)

def test(space):
    if test_numpy:
        global printMsg
        printMsg = True

        numpySpace  = create.space(space, grid, dimRange=1, order=1, storage='numpy')
        numpy_h     = create.function("discrete", numpySpace, name="numpy")
        numpy_dofs  = numpy_h.as_numpy

        def preconditioner(u , v):
            global printMsg
            if printMsg:
                print(f"Python fct preconditioner: u = {u.name} and v = {v.name}")
            printMsg = False

            # do nothing preconditioner
            v.assign(u)
            return

        code = """
        #include <iostream>
        template< class DF >
        void preconditioner( const DF &u, DF &v ) {
          v.assign(u);
        }
        """
        preconditioner = algorithm.load('preconditioner', io.StringIO(code), numpy_h,numpy_h)

        solverParameters = { "newton.linear.preconditioning.method": preconditioner }
        numpyScheme = create.scheme("galerkin", model, numpySpace, parameters=solverParameters)

        numpyScheme.solve(numpy_h)
        start= time.time()
        for i in range(testLoop):
            linOp   = linearOperator(numpyScheme)
            numpyScheme.jacobian(numpy_h,linOp)
            numpy_mat = linOp.as_numpy
        end = time.time()
        # print( "numpy:", (end-start)/testLoop, flush=True )
        # sys.stdout.flush()
        numpy_coo   = numpy_mat.tocoo()
        # for i,j,v in zip(numpy_coo.row,numpy_coo.col,numpy_coo.data):
        #     print(i,j,v)
        # print("****************************",flush=True)

    if test_istl:
        istlSpace  = create.space(space, grid, dimRange=1, order=1, storage='istl')
        istlScheme = create.scheme("galerkin", model, istlSpace)
        istl_h     = create.function("discrete", istlSpace, name="istl")
        istl_dofs  = istl_h.as_istl
        # istlScheme.solve(target = istl_h)
        start= time.time()
        for i in range(testLoop):
            linOp = linearOperator(istlScheme)
            istlScheme.jacobian(istl_h,linOp)
            istl_mat = linOp.as_istl
        end= time.time()
        # print( "istl:", (end-start)/testLoop, flush=True )
        # sys.stdout.flush()
        # there is no way yet to go from istl to scipy - would be nice to have
        # istl_coo = istl_mat.tocoo()
        # for i,j,v in zip(eigen_coo.row,eigen_coo.col,eigen_coo.data):
        #     print(i,j,v)
        # print("****************************",flush=True)

    if test_petsc:
        import petsc4py
        from petsc4py import PETSc
        petsc4py.init(sys.argv)
        petscSpace  = create.space(space, grid, dimRange=1, order=1, storage='petsc')
        petscScheme = create.scheme("galerkin", model, petscSpace)
        petsc_h     = create.function("discrete", petscSpace, name="petsc")
        petsc_dofs  = petsc_h.as_petsc
        # petscScheme.solve(target = petsc_h)
        linOp = linearOperator(petscScheme)
        petscScheme.jacobian(petsc_h,linOp)
        petsc_mat = linOp.as_petsc
        rptr, cind, vals = petsc_mat.getValuesCSR()
        petsc_coo = scipy.sparse.csr_matrix((vals,cind,rptr)).tocoo()
        start= time.time()
        for i in range(testLoop):
            linOp = linearOperator(petscScheme)
            petscScheme.jacobian(petsc_h,linOp)
            petsc_mat = linOp.as_petsc
        end= time.time()
        # print( "petsc:", (end-start)/testLoop, flush=True )
        # sys.stdout.flush()

        ksp = PETSc.KSP()
        ksp.create(PETSc.COMM_WORLD)
        ksp.setType("cg")
        # ksp.getPC().setType("icc")
        petsc_h.clear()
        res = petsc_h.copy()
        petscScheme(petsc_h, res)
        petscScheme.jacobian(petsc_h,linOp)
        petsc_mat = linOp.as_petsc
        ksp.setOperators(petsc_mat,petsc_mat)
        ksp.setFromOptions()
        ksp.solve(res.as_petsc, petsc_h.as_petsc)


        # print("****************************",flush=True)
        # print(petsc_mat.size, petsc_mat.getSize(), petsc_mat.getSizes())
        # print(petsc_mat.getType())
        # print(type(petsc_mat))
        # print(petscSpace.size)
        # print(petsc_mat.assembled)
        # rptr, cind, vals = petsc_mat.getValuesCSR()
        # petsc_coo = scipy.sparse.csr_matrix((vals,cind,rptr),shape=(100,100)).tocoo()
        # for i,j,v in zip(petsc_coo.row,petsc_coo.col,petsc_coo.data):
        #     print(i,j,v)

    if test_istl:
        try: # istl_coo does not exist
            assert (istl_coo.row == numpy_coo.row).all()
            assert (istl_coo.col == numpy_coo.col).all()
            assert np.allclose( istl_coo.data, numpy_coo.data )
        except:
            # print("issue between istl and numpy matrices")
            pass
    if test_petsc:
        try:
            assert (petsc_coo.row == numpy_coo.row).all()
            assert (petsc_coo.col == numpy_coo.col).all()
            assert np.allclose( petsc_coo.data, numpy_coo.data )
        except:
            # print("issue between petsc and numpy matrices")
            pass

test("lagrange")
# test("dglagrange")
# test("dglegendre") # only works with cube grids
