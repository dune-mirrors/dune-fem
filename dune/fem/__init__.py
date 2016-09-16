from __future__ import print_function
from mpi4py import MPI
from ..femmpi import comm
from . import grid
import dune.models.gridfunction
from ..common import reader
from .create import *

def string2dgf(dgf):
    return (reader.dgfString,"DGF\n" + dgf)
def cartesianDomain(lower,upper,division):
    ret = "INTERVAL\n"
    for x in lower:
        ret = ret + str(x) + " "
    ret = ret + "\n"
    for x in upper:
        ret = ret + str(x) + " "
    ret = ret + "\n"
    for x in division:
        ret = ret + str(x) + " "
    ret = ret + "\n#\n"
    return string2dgf(ret)
