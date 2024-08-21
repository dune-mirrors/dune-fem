import time
import numpy as np
import ufl
import dune.ufl
import dune.grid
import dune.fem

# %%
dimRange   = 1
dt         = 0.1
spiral_a   = 0.75
spiral_b   = 0.02
spiral_eps = 0.02
spiral_D   = 1./100
def spiral_h(u,v): return u - v

x = ufl.SpatialCoordinate(dune.ufl.domain(2))
initial_u = ufl.conditional(x[1]>1.25,1,0)
initial_v = ufl.conditional(x[0]<1.25,0.5,0)

solverParameters = {
        "nonlinear.tolerance": 1e-3,
        "nonlinear.verbose": False,
        "linear.tolerance": 1e-5,
        "linear.preconditioning.method": "none",
        "linear.verbose": False}

gridView = dune.grid.structuredGrid([0,0],[2.5,2.5],[100,100])

def simulate(spcStorage,solver):
    space = dune.fem.space.lagrange( gridView, dimRange=dimRange, order=1,
                                     storage=spcStorage )

    uh   = space.interpolate( initial_u, name="u" )
    uh_n = uh.copy()
    vh   = space.interpolate( initial_v, name="v" )
    vh_n = vh.copy()

    u   = ufl.TrialFunction(space)
    phi = ufl.TestFunction(space)

    a_ex = ufl.inner(uh_n, phi) * ufl.dx
    a_im = (dt * spiral_D * ufl.inner(ufl.grad(u), ufl.grad(phi)) +
            ufl.inner(u,phi)) * ufl.dx

    ustar = (vh_n[0]+spiral_b)/spiral_a
    a_ex += ufl.conditional(uh_n[0]<ustar, dt/spiral_eps* u[0]*(1-uh_n[0])*(uh_n[0]-ustar),
                                           dt/spiral_eps*uh_n[0]*(1-u[0]) *(uh_n[0]-ustar) ) * phi[0] * ufl.dx

    equation   = a_im == a_ex
    ode_update = ufl.as_vector([ vh_n[0] + dt*spiral_h(uh_n[0], vh_n[0]) ])
    # make sure the gf has been compiled before the time loop
    vh_n.interpolate( ode_update )

    scheme = dune.fem.scheme.galerkin( equation, space, parameters=solverParameters,
                                       solver=solver)

    # should be no compilation left after this point...
    print(f"starting simulation with {spcStorage},{solver}",flush=True)
    t = 0.
    linIter = 0
    nlinIter = 0

    start = time.time()
    while t < 10:
        uh_n.assign(uh)
        vh_n.assign(vh)
        info = scheme.solve(target=uh)
        linIter += info["linear_iterations"]
        nlinIter += info["iterations"]
        vh.interpolate( ode_update )
        t += dt
    used = time.time() - start
    print(f"    time used: {used}",flush=True)
    print(f"    iterations: {nlinIter}, {linIter}",flush=True)

    # uh.plot()
    return uh,vh, used

uref,vref,_ = simulate(spcStorage="numpy",solver="cg")

uh,vh,tnp = simulate(spcStorage="numpy",solver=("petsc","cg"))
e_u = uref.as_numpy - uh.as_numpy
e_v = vref.as_numpy - vh.as_numpy
print("    test:",np.dot(e_u,e_u),np.dot(e_v,e_v),flush=True)

uh,vh,tpp = simulate(spcStorage="petsc",solver="cg")
e_u = uref.as_numpy - uh.as_petsc.getArray()
e_v = vref.as_numpy - vh.as_petsc.getArray()
print("    test:",np.dot(e_u,e_u),np.dot(e_v,e_v),flush=True)

print("FINAL:", abs(tnp-tpp)/tnp,flush=True)
