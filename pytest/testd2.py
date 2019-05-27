from __future__ import print_function

import math
from ufl import *
import ufl

from dune.grid import cartesianDomain
from dune.fem import parameter
from dune.ufl import Space
from dune.fem.function import integrate

import dune.create as create

parameter.append({"fem.verboserank": 0})
newtonParameter = {"tolerance": 1e-7, "verbose": "false",
                   "linear.absolutetol": 1e-8, "linear.reductiontol": 1e-8,
                   "linear.preconditioning.method": "jacobi",
                   "linear.preconditioning.iterations": 1, "linear.preconditioning.relaxation": 1.2,
                   "linear.verbose": "false"}

grid = create.grid("ALUConform", cartesianDomain([0,0],[1,1],[10,10]), dimgrid=2)
spc = create.space("lagrange", grid, dimRange=1, order=2, storage="istl")

uflSpace = Space((grid.dimGrid, grid.dimWorld), 1)
u = TrialFunction(uflSpace)
v = TestFunction(uflSpace)
x = SpatialCoordinate(uflSpace.cell())
n = FacetNormal(uflSpace.cell())
mu = 1.
# he = MaxFacetEdgeLength(uflSpace.cell())('+') # this is wrong
# hT = MaxCellEdgeLength(uflSpace.cell())
hT = CellVolume(uflSpace.cell())
hF = FacetArea(uflSpace.cell())
# he = FacetArea(uflSpace.cell()) / Min( avg('+'), avg('-') )
heInv = hF / avg( hT )
exact = as_vector( [cos(pi*x[0])*cos(pi*x[1])] )

#########

a  = inner(grad(u-exact), grad(v)) * dx
a += mu*hT * div(grad(u[0]-exact[0])) * div(grad(v[0])) * dx
s  = mu/heInv * inner( jump(grad(u[0])), jump(grad(v[0])) ) * dS
s += mu/hF * inner( u-exact, v ) * ds
model  = create.model("integrands", grid, a+s == 0)
scheme = create.scheme("galerkin", model, spc, solver="cg",
            parameters={"newton." + k: v for k, v in newtonParameter.items()})
solA = spc.interpolate([0],name="solA")
scheme.solve(solA)

########

a  = div(grad(u[0]-exact[0])) * div(grad(v[0])) * dx
s  = mu*heInv * inner( jump(grad(u[0])), jump(grad(v[0])) ) * dS
s += mu/hF**3 * inner( u-exact, v ) * ds
model  = create.model("integrands", grid, a+s == 0)
scheme = create.scheme("galerkin", model, spc, solver="cg",
                parameters={"newton." + k: v for k, v in newtonParameter.items()})
solB = spc.interpolate([0],name="solB")
scheme.solve(solB)

errA_sol = math.sqrt( integrate(grid, (solA-exact)**2, 5) )
errB_sol = math.sqrt( integrate(grid, (solB-exact)**2, 5) )
errA_B   = math.sqrt( integrate(grid, (solA-solB)**2, 5) )

print( errA_sol, errB_sol, errA_B ) # 0.0004520603651576     0.013241522498765897 0.012944687615068362

assert abs(errA_sol-0.00045)/0.00045 < 0.1
assert abs(errB_sol-0.013)/0.013 < 0.1
assert abs(errA_B-0.013)/0.013 < 0.1
