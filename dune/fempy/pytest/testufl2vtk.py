import ufl
x = ufl.SpatialCoordinate(ufl.triangle)
u=ufl.sin(x[0]*x[1])
import dune.grid
grid=dune.grid.structuredGrid([0,0],[1,1],[10,10])
import dune.create as create
test = create.function("ufl", grid, "cutoff", 5,
       ufl.as_vector([u]))
grid.writeVTK("test", pointdata={test.name:test, "ufl":u},
                      pointvector={"grad":ufl.grad(u)} )
