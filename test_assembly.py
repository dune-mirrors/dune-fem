import time,sys
import dune
import dune.create as create

try:
    from dune.fem.operator import linear as linearOperator
    def jacobian(scheme,uh,A):
        scheme.jacobian(uh,A)
except: # for 2.6
    def linearOperator(scheme):
        u_h    = create.function("discrete", scheme.space, name="tmp")
        u_h.clear()
        return scheme.assemble(u_h)
    def jacobian(scheme,uh,A):
        A = scheme.assemble(uh)


from dune.ufl import Space
from ufl import Identity, TestFunction, TrialFunction, SpatialCoordinate, ds, dx, inner, grad, div

test_fem   = True
test_istl  = True
test_petsc = True

test_scalar = True
test_vector = True
test_21 = True
test_12 = True
testLoop = 10



grid = create.grid("ALUConform", dune.grid.cartesianDomain([0, 0], [1, 1], [189, 189]), dimgrid=2)

def test(model,spaceName,dimD,dimR,storage):
    print("########################################")
    print("#### ",spaceName,storage,dimD,dimR,flush=True)
    spaceD  = create.space(spaceName, grid, dimRange=dimD, order=1, storage=storage)
    spaceR  = create.space(spaceName, grid, dimRange=dimR, order=1, storage=storage)
    scheme = create.operator("galerkin", model, spaceD, spaceR)
    uD     = create.function("discrete", spaceD, name=storage)
    uD.clear()
    start = time.clock()
    for i in range(testLoop):
        A = linearOperator(scheme)
    end = time.clock()
    print( "setup+assembly:",(end-start)/testLoop, flush=True )
    start = time.clock()
    for i in range(testLoop):
        jacobian(scheme,uD,A)
    end = time.clock()
    print( "assembly only: ",(end-start)/testLoop, flush=True )
    sys.stdout.flush()

    try:
        import petsc4py
        from petsc4py import PETSc
        mat = A.as_petsc
        print(mat.getInfo(), flush=True)
    except:
        pass
    print("########################################")

if test_scalar:
    dimRange = 1
    uflSpace = Space((2,2),dimRange)
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())

    rhs = (x[0] + x[1]) * v[0]
    a = inner(grad(u), grad(v)) * dx
    b = rhs * dx
    model = create.model("integrands", grid, a==b)

    space = "lagrange"
    if test_fem:
        test(model,space,dimRange,dimRange,"fem")
    if test_istl:
        test(model,space,dimRange,dimRange,"istl")
    if test_petsc:
        test(model,space,dimRange,dimRange,"petsc")

if test_vector:
    dimRange = 2
    uflSpace = Space((2,2),dimRange)
    u = TrialFunction(uflSpace)
    v = TestFunction(uflSpace)
    x = SpatialCoordinate(uflSpace.cell())

    rhs = (x[0] + x[1]) * v[0]
    a = inner(grad(u), grad(v)) * dx
    b = rhs * dx
    model = create.model("integrands", grid, a==b)

    space = "lagrange"
    if test_fem:
        test(model,space,dimRange,dimRange,"fem")
    if test_istl:
        test(model,space,dimRange,dimRange,"istl")
    if test_petsc:
        test(model,space,dimRange,dimRange,"petsc")

#############################
#############################
#############################

uflSpace1 = Space((2,2),1)
u1 = TrialFunction(uflSpace1)
v1 = TestFunction(uflSpace1)
uflSpace2 = Space((2,2),2)
u2 = TrialFunction(uflSpace2)
v2 = TestFunction(uflSpace2)

if test_21:
    a = div(u2)*v1[0]*dx
    model = create.model("integrands", grid, a)

    space = "lagrange"
    if test_fem:
        test(model,space,2,1,"fem")
    if test_istl:
        test(model,space,2,1,"istl")
    if test_petsc:
        test(model,space,2,1,"petsc")

if test_12:
    a = -inner( u1[0]*Identity(2), grad(v2) ) * dx
    model = create.model("integrands", grid, a)

    space = "lagrange"
    if test_fem:
        test(model,space,1,2,"fem")
    if test_istl:
        test(model,space,1,2,"istl")
    if test_petsc:
        test(model,space,1,2,"petsc")
