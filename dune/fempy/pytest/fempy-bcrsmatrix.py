import dune.grid
import dune.create as create
from dune.istl import BlockVector
#create a blockvector which composes of 1 Fieldvector of size 2
vec = BlockVector(2)(10)

# print(vec, len(vec), vec[0], vec[0][0])
#set first element in first field vector equal to one
vec[0][0] = 2
# print(vec, len(vec), vec[0], vec[0][0])

#change element to 1
vec[0][0] = 1

assert(vec[0][0] == 1),"dune.istl blockvector not changing value correctly "
#now test with with other fieldvector created in dune=fempy

view = create.grid("yasp", dune.grid.cartesianDomain([0,0],[1,1],[1,1]),dimgrid=2)

spc = create.space("Lagrange", view, dimRange=2, order=3, storage='istl')

#initialise all to 2
uh = spc.interpolate([2,3], name="solution")

df = uh.as_istl
# print(df, len(df), df[0], df[0][0])
df[0][0] = 1
# print(df, len(df), df[0], df[0][0])

# print(uh.dofVector[0][0])

#reset equalt to one againt
assert(df[0][0] == 1),"blockvector from dune-fempy not changing correctly"
