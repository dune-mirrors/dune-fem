from __future__ import print_function

import dune.common

from ._fem import *
from ._adaptation import adapt, loadBalance

from . import view as view
from . import space as space
from . import discretefunction as discretefunction
from . import scheme as scheme
from . import function as function
from . import model as model

registry = {}

registry["view"] = {
         "adaptive"   : view.adaptiveLeafGridView,
         "filtered"   : view.filteredGridView,
         "geometry"   : view.geometryGridView
     }
registry["space"] = {
         "Lagrange"   : space.lagrange,
         "DGONB"      : space.dgonb,
         "P1Bubble"   : space.p1Bubble
     }
registry["discretefunction"] = {
         "adaptive" : discretefunction.adaptive,
         "fem"      : discretefunction.adaptive,
         "istl"     : discretefunction.istl,
         "eigen"    : discretefunction.eigen
     }
registry["scheme"] = {
         "h1"      : scheme.h1,
         "dg"      : scheme.dg,
         "mvdg"    : scheme.nvdg,
         "stokes"  : scheme.stokes,
         "burgers" : scheme.burgers
     }
registry["function"] = {
         "global"     : function.globalFunction,
         "local"      : function.localFunction,
         "cpp"        : function.cppFunction,
         "ufl"        : function.uflFunction,
         "numpy"      : function.numpyFunction,
         "levels"     : function.levelFunction,
         "partitions" : function.partitionFunction,
         "discrete"   : function.discreteFunction
     }
registry["model"] = {
         "elliptic"   : model.elliptic
     }
