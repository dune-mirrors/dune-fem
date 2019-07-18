from __future__ import print_function, unicode_literals
import time, math, numpy, ufl
import dune, dune.create
from dune.generator import algorithm

from dune.grid import cartesianDomain, yaspGrid
domain = cartesianDomain([0, 0], [1, 0.25], [12, 3])
yaspView = yaspGrid(domain)

x = ufl.SpatialCoordinate( ufl.triangle )
function = ufl.as_vector([ ufl.cos(2*ufl.pi/(0.3+x[0]*x[1])) ])
# function = dune.create.function("global",gridview=yaspView, name="gl",order=5,
#            value=lambda x: [math.cos(2*math.pi/(0.3+x[0]*x[1]))] )
space = dune.create.space("lagrange", yaspView, order=1, dimRange=1)
uh = space.interpolate( function, name="uh" )
error = dune.create.function("ufl",gridView=yaspView,name="error",order=5,ufl=uh-function,
        virtualize=True)

rules = dune.geometry.quadratureRules(5)

if False: # dune.fem type grid functions don't yet have vectorization support
    start = time.time()
    l2norm2 = 0
    for e in yaspView.elements:
        hatxs, hatws = rules(e.type).get()
        weights = hatws * e.geometry.integrationElement(hatxs)
        l2norm2 += numpy.sum(error(e, hatxs)**2 * weights, axis=-1)
    print("Python:",math.sqrt(l2norm2),flush=True)
    print("time used:", round(time.time()-start,2),flush=True)

if True:
    #algo = algorithm.load('l2norm2', 'test_quad.hh', yaspView, rules, error)
    algo = algorithm.load('l2norm2FemQuad', 'test_quad.hpp', yaspView, rules, error)
    start = time.time()
    l2norm2 = algo(yaspView,rules,error)
    print("C++:",math.sqrt(l2norm2),flush=True)
    print("time used:", round(time.time()-start,2),flush=True)

try:
    import dune.geometry.quadpy as quadpy
    rules = quadpy.rules({dune.geometry.quadrilateral: ("C2 7-2","Stroud")})

    if False:
        start = time.time()
        l2norm2 = 0
        for e in yaspView.elements:
            hatxs, hatws = rules(e.type).get()
            weights = hatws * e.geometry.integrationElement(hatxs)
            l2norm2 += numpy.sum(error(e, hatxs)**2 * weights, axis=-1)
        print("Python:",math.sqrt(l2norm2),flush=True)
        print("time used:", round(time.time()-start,2),flush=True)

    if True:
        start = time.time()
        l2norm2 = algo(yaspView,rules,error)
        print("C++:",math.sqrt(l2norm2),flush=True)
        print("time used:", round(time.time()-start,2),flush=True)
except ImportError:
    pass
