#!/usr/bin/env python3

from dune.fem.space import finiteVolume as solSpace
from dune.grid import cartesianDomain
from dune.grid import yaspGrid as leafGridView

domain = cartesianDomain([0, 0], [1, 1], [10, 10])
grid = leafGridView(domain)
grid.hierarchicalGrid.globalRefine(2)

views = [grid.hierarchicalGrid.levelView(i) for i in range(3)]
spaces = [solSpace(view) for view in views]

spc0 = solSpace(views[0])
spc1 = solSpace(views[1])

u = spc0.DiscreteFunction(spc0,'u')
assert u.space.size == u.as_numpy.size, f"size of u {u.space.size} and dofvec {u.as_numpy.size} do not match"
assert u.space.size == spc0.size, f"size of u {u.space.size} and space {spc0.size} do not match"

# works and has the correct sizes
for space in spaces:
    assert ( space.size == space.gridView.size(0) )

# does not work, sizes are mixed up
for space in spaces:
    u = space.function(name='u')
    assert u.as_numpy.size == space.size, f"Size of space {space.size} and function {u.as_numpy.size} differs"
    assert u.space.size == space.size, f"Size of space {space.size} and function {u.space.size} differs"
