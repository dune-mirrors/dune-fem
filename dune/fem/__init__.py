from __future__ import print_function
from mpi4py import MPI
from .grid import leafGrid
from .space import create as space

def interpolate(grid, func, **kwargs):
    storage=kwargs.pop('a',"Adaptive")
    # try:
    #     storage=kwargs['storage']
    # except:
    #     storage="Adaptive"

    # spaceName = kwargs['space']
    spaceName = kwargs.pop('space')
    mySpace=space(spaceName,grid,**kwargs)
    return mySpace.interpolate(func,**kwargs)
