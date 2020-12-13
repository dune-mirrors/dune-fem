"""
FEniCS tutorial demo program: Linear elastic problem.

  -div(sigma(u)) = f

The model is used to simulate an elastic beam clamped at
its left end and deformed under its own weight.
"""

from __future__ import print_function
from dune.fenics import *

# Scaled variables
L = 1; W = 0.2
mu = 1
rho = 1
delta = W/L
gamma = 0.4*delta**2
beta = 1.25
lambda_ = beta
g = gamma

# Create mesh and define function space
mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 10, 3, 3)
V = VectorFunctionSpace(mesh, 'P', 1)

# Define boundary condition
tol = 1E-14

# def clamped_boundary(x, on_boundary):
#     return on_boundary and x[0] < tol

# bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)
# bc = DirichletBC(V, Constant((0,0,0)), 1)
bc = DirichletBC(V, as_vector([0,0,0]), 1)

# Define strain and stress

def epsilon(u):
    return 0.5*(nabla_grad(u) + nabla_grad(u).T)
    #return sym(nabla_grad(u))

def sigma(u):
    return lambda_*nabla_div(u)*Identity(d) + 2*mu*epsilon(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0, 0, -rho*g))
T = Constant((0, 0, 0))
a = inner(sigma(u), epsilon(v))*dx
# L = dot(f, v)*dx + dot(T, v)*ds
L = dot(as_vector([0,0,-rho*g]),v)*dx

# Compute solution
u = Function(V)
solve(a == L, u, bc)

# Plot solution (can't yet do 3d - what can fenics do?)
# plot(u) # , title='Displacement', mode='displacement')

# Plot stress
s = sigma(u) - (1./3)*tr(sigma(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
V = FunctionSpace(mesh, 'P', 1)
# von_Mises = project(von_Mises, V)
von_Mises = interpolate(von_Mises, V)
plot(von_Mises, title='Stress intensity')

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
# u_magnitude = project(u_magnitude, V)
u_magnitude = interpolate(u_magnitude, V)
# plot(u_magnitude, 'Displacement magnitude')
print('min/max u:',
      u_magnitude.pointData().min(),
      u_magnitude.pointData().max())

# Save solution to file in VTK format
# File('elasticity/displacement.pvd') << u
# File('elasticity/von_mises.pvd') << von_Mises
# File('elasticity/magnitude.pvd') << u_magnitude
mesh.writeVTK("fenics-elasticity", pointdata=[u,von_Mises,u_magnitude])

# Hold plot
# interactive()
