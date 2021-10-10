from __future__ import print_function, division
import math
import numpy

from dune.grid import structuredGrid
import dune.create as create

grid = structuredGrid([0,0],[1,1],[10,10])

fvspc = create.space("finitevolume", grid, dimRange=2, storage="numpy")
estimate = fvspc.interpolate([1,1], name="estimate")

lf = estimate.localFunction(grid.elements.__next__())
y = lf.evaluate([0,0])
assert y == (1,1)
# print(y)

lf = estimate.setLocalContribution()
lf.bind(grid.elements.__next__())
y1,y2 = lf[0],lf[1]
assert y1 == 0 and y2 == 0
# print(y1,y2)
lf[0] = 10
y1,y2 = lf[0],lf[1]
assert y1 == 10 and y2 == 0
# print(y1,y2)
lf.unbind()

lf = estimate.localContribution("set")
lf.bind(grid.elements.__next__())
y1,y2 = lf[0],lf[1]
assert y1 == 0 and y2 == 0
# print(y1,y2)
lf[0] = 10
lf[1] = -10
y1,y2 = lf[0],lf[1]
assert y1 == 10 and y2 == -10
# print(y1,y2)
lf.unbind()

lf = estimate.localFunction(grid.elements.__next__())
y = lf.evaluate([0,0])
assert y == (10,-10)
# print(y)
