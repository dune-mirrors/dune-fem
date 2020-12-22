import math
from ufl import SpatialCoordinate, sin, pi, dot, as_vector

from dune.grid import structuredGrid, gridFunction
from dune.fem.space import lagrange

grid = structuredGrid([0,0],[1,1],[5,5])
space = lagrange(grid, order=3, dimRange=2)
x = SpatialCoordinate(space)
df = space.interpolate([ sin(2**i*pi*x[0]*x[1]/(0.25**(2*i)+dot(x,x)))
                         for i in range(1,3)],
                       name="femgf")

@gridFunction(grid)
def f(x):
    return [ math.sin(2**i*math.pi*x[0]*x[1]/(0.25**(2*i)+x[0]*x[0]+x[1]*x[1]))
             for i in range(1,3) ]

from dune.vtk import vtkWriter
writer = vtkWriter( grid, "testFem",
                    pointScalar = {  "f"        : f,
                                     "df"       : df,
                                    ("f0",(1,)) : f,
                                     "df0"      : df[1],
                                     "sin"      : [ sin(2**i*pi*x[0]*x[1]/(0.25**(2*i)+dot(x,x)))
                                                    for i in range(1,3) ]
                                  },
                    pointData   = [df]
                  )
writer = vtkWriter( grid, "testFem_lag3",
                    version="lagrange", order=3,
                    pointData   = { "dfV":    df,
                                    "df3D":   as_vector([df[0],df[1],abs(dot(df,x))]),
                                    "df1V":   df[1],
                                    "tensor": as_vector([df[0],df[1],df[0],df[1]])
                                  },
                    pointScalar = {  "f"        : f,
                                     "df"       : df,
                                    ("f1",(1,)) : f,
                                     "df1"      : df[1],
                                     "sin"      : [ sin(2**i*pi*x[0]*x[1]/(0.25**(2*i)+dot(x,x)))
                                                    for i in range(1,3) ]
                                  },
                  )
