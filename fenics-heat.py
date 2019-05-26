"""
FEniCS tutorial demo program: Diffusion of a Gaussian hill.

  u'= Laplace(u) + f  in a square domain
  u = u_D             on the boundary
  u = u_0             at t = 0

  u_D = f = 0

The initial condition u_0 is chosen as a Gaussian hill.
"""
from __future__ import print_function
from dune.fenics import *
import time
import dune.fem
dune.fem.parameter.append({"fem.verboserank": 0})


T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size

# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
# def boundary(x, on_boundary):
#     return on_boundary
g_D = Expression("t/endTime*sin(pi/2*x[0]*x[1])",mesh=mesh,degree=5,t=0,endTime=T)
print(g_D.endTime)
bc = [DirichletBC(V, g_D, i) for i in range(1,5)] # boundary ids are 1,..,4 or yasp
# bc = [DirichletBC(V, Constant(0), i) for i in range(1,5)] # boundary ids are 1,..,4 or yasp

# Define initial value
u_0 = Expression('exp(-a*pow(x[0], 2) - a*pow(x[1], 2)) + g_D', mesh=mesh,
                 degree=2, a=5, g_D=g_D)
u_n = interpolate(u_0, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Create VTK file for saving solution
# vtkfile = File('heat_gaussian/solution.pvd')

# Time-stepping
u = Function(V)
t = 0
vtkfile = mesh.sequencedVTK("fenics-heat", pointdata=[u])
for n in range(num_steps):

    # Update current time
    t += dt
    g_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Save to file and plot solution
    # vtkfile << (u, t)
    if n%10 == 0:
        vtkfile()
        plot(u)

    # Update previous solution
    u_n.assign(u)

# Hold plot
# interactive()
