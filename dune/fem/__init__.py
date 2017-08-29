from __future__ import print_function

import dune.common
from ._fem import *
from ._adaptation import adapt, loadBalance

from . import view as view
from . import space as space
from . import discretefunction as discretefunction
from . import operator as operator
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
         "lagrange"      : space.lagrange,
         "dgonb"         : space.dgonb,
         "dglegendre"    : space.dglegendre,
         "dglagrange"    : space.dglagrange,
         "fv"            : space.finiteVolume,
         "p1bubble"      : space.p1Bubble,
         "combined"      : space.combined,
         "tuple"         : space.combined,
         "bdm"           : space.bdm,
         "rannacherTurek": space.rannacherTurek
     }
registry["discretefunction"] = {
         "adaptive" : discretefunction.adaptive,
         "fem"      : discretefunction.adaptive,
         "istl"     : discretefunction.istl,
         "eigen"    : discretefunction.eigen,
         "petsc"    : discretefunction.petsc
     }
registry["operator"] = {
        "galerkin"   : operator.galerkin
    }
registry["solver"] = {
         "fem"         : discretefunction.femsolver,
         "pardg"       : discretefunction.pardgsolver,
         "femoem"      : discretefunction.oemfemsolver,
         "istl"        : discretefunction.istlsolver,
         "suitesparse" : discretefunction.suitesparsesolver,
         "eigen"       : discretefunction.eigensolver,
         "viennacl"    : discretefunction.viennaclsolver,
         "petsc"       : discretefunction.petscsolver
     }
registry["scheme"] = {
         "h1"         : scheme.h1,
         "h1galerkin" : scheme.h1Galerkin,
         "dg"         : scheme.dg,
         "dggalerkin" : scheme.dgGalerkin,
         "galerkin"   : scheme.galerkin,
         "linearized" : scheme.linearized,
         "mvdg"       : scheme.nvdg,
         "stokes"     : scheme.stokes,
         "burgers"    : scheme.burgers
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
         "elliptic"   : model.elliptic,
         "integrands" : model.integrands
     }

from ufl import dx
from dune.ufl import GridFunction
def evaluate(expression,x,**kwargs):
    # don't know how to get the coefficients form an expression
    coefficients = set((expression**2*dx).coefficients())
    for coefficient in coefficients:
        if coefficient.is_cellwise_constant():
            try:
                name  = coefficient.str()
                value = kwargs.get(name,None)
                if value:
                    pass
                    # replace the coefficient with the value
            except AttributeError:
                pass
        # else:
        #     if isinstance(coefficient,GridFunction):
        #         coefficient.bind(x)
    val = expression(x)
    for coefficient in coefficients:
        if coefficient.is_cellwise_constant():
            pass
        else:
            if isinstance(coefficient,GridCoefficient):
                coefficient.unbind()
    return val
