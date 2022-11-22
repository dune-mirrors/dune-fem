from __future__ import print_function
import atexit

import dune.common
from ._fem import *
from ._adaptation import adapt, loadBalance, mark, markNeighbors, globalRefine, doerflerMark
from ._spaceadaptation import spaceAdapt

from . import view as view
from . import space as space
from . import discretefunction as discretefunction
from . import operator as operator
from . import scheme as scheme
from . import function as function
from . import model as model
from . import plotting

# finalization of fem module (i.e. calling PETSc finalize etc)
# atexit.register( _fem.__finalizeFemModule__ )

registry = {}

registry["view"] = {
         "adaptive"   : view.adaptiveLeafGridView,
         "filtered"   : view.filteredGridView,
         "geometry"   : view.geometryGridView
     }
registry["space"] = {
         "lagrange"      : space.lagrange,
         "lagrangehp"    : space.lagrangehp,
         "dgonb"         : space.dgonb,
         "dgonbhp"       : space.dgonbhp,
         "dglegendre"    : space.dglegendre,
         "dglegendrehp"  : space.dglegendrehp,
         "dglagrange"    : space.dglagrange,
         "finitevolume"  : space.finiteVolume,
         "p1bubble"      : space.p1Bubble,
         "combined"      : space.combined,
         "product"       : space.product,
         "bdm"           : space.bdm,
         "raviartthomas" : space.raviartThomas,
         "rannacherTurek": space.rannacherTurek
     }
registry["discretefunction"] = {
         "numpy"      : discretefunction.numpy,
         "adaptive"   : discretefunction.adaptive,
         "fem"        : discretefunction.fem,
         "istl"       : discretefunction.istl,
         "eigen"      : discretefunction.eigen,
         "petsc"      : discretefunction.petsc,
         "petscadapt" : discretefunction.petscadapt
     }
registry["operator"] = {
        "galerkin"   : operator.galerkin,
        "h1"         : operator.h1
    }
registry["solver"] = {
         "fem"         : discretefunction.femsolver,
         "istl"        : discretefunction.istlsolver,
         "suitesparse" : discretefunction.suitesparsesolver,
         "eigen"       : discretefunction.eigensolver,
         "viennacl"    : discretefunction.viennaclsolver,
         "petsc"       : discretefunction.petscsolver,
         "amgx"        : discretefunction.amgxsolver
     }
registry["scheme"] = {
         "h1"         : scheme.h1,
         "h1galerkin" : scheme.h1Galerkin,
         "dg"         : scheme.dg,
         "dggalerkin" : scheme.dgGalerkin,
         "galerkin"   : scheme.galerkin,
         "linearized" : scheme.linearized
     }
registry["function"] = {
         "global"     : function.globalFunction,
         "local"      : function.localFunction,
         "cpp"        : function.cppFunction,
         "ufl"        : function.uflFunction,
         "levels"     : function.levelFunction,
         "partitions" : function.partitionFunction,
         "discrete"   : function.discreteFunction
     }
registry["model"] = {
         "elliptic"   : model.elliptic,
         "integrands" : model.integrands
     }

#------------------------------------------------------
# Setting verbosity level for dune-fem
def setVerbosity( level=1, rank=0 ):
    """
    Verbosity levels:
      0: no output
      1: solver statistics (default)
      2: extended solver statistics
      3: parameter and other output
      4: diagnostics output
      5: debug output

    Set verbose rank (-1,0,...,MPI_size-1) to the rank that should print output.
    """
    parameter._setVerbosity( level, rank )
#------------------------------------------------------

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

from dune.grid.grid_generator import _writeVTKDispatcher
def vtkDispatchUFL(grid,f):
    from dune.fem.function._functions import uflFunction
    order = 5 # needs to be derived from f
    # gf = uflFunction(grid, "tmp", order, f)
    try:
        gf = uflFunction(grid, "tmp", order, f).gf
    except AttributeError:
        gf = None
    return gf
_writeVTKDispatcher.append(vtkDispatchUFL)

import pickle as _pickle
import pickletools
import importlib, sys
def dump(objects, f, protocol=None):
    if not (isinstance(objects,list) or isinstance(objects,tuple)):
        raise TypeError("only pass in tuples/lists to pickle")
    gv = None
    obj = []
    for o in objects:
        if hasattr(o,"gridView"):
            if gv and not (gv==o.gridView):
                raise TypeError("all gridviews must be the same")
            gv = o.gridView
        if hasattr(o,"__impl__"):
            obj.append((True,o.__impl__))
        else:
            obj.append((False,o))
    objDump = _pickle.dumps(obj,protocol)
    mods = []
    for opcode,arg,pos in pickletools.genops(objDump):
        try:
            if "dune.generated" in arg:
                mods += [arg]
        except:
            pass
    modsCode = []
    for m in mods:
        fname = importlib.import_module(m).__file__
        fname = fname.replace(".so",".cc")
        with open(fname,"r") as fmod:
            modsCode += [[m, fmod.read()]]
    _pickle.dump(modsCode,f,protocol)
    _pickle.dump(obj,f,protocol)

def load(f):
    # make sure dune.generated in dune-py is added to dune package path
    from dune.generator import builder
    builder.initialize()

    """
    mods = []
    for opcode,arg,pos in pickletools.genops(f):
        try:
            if "dune.generated" in arg:
                mods += [arg]
        except:
            pass
    while len(mods)>0:
        cont = False
        for m in mods:
            module = sys.modules.get(m)
            try:
                importlib.import_module(m)
                print("loading",m)
                mods.remove(m)
                cont = True
            except ImportError as ex:
                print(ex)
                print("skipping",m)
                pass
        assert cont, "no module could be loaded"
    f.seek(0)
    """
    mods = _pickle.load(f)
    while len(mods)>0:
        cont = False
        for m in mods:
            try:
                module = builder.load(m[0].split(".")[-1], m[1], None)
                # print("loading",m[0])
                mods.remove(m)
                cont = True
            except ImportError as ex:
                # print(ex)
                # print("skipping",m[0])
                pass
        assert cont, "no module could be loaded"

    objects = []
    obj = _pickle.load(f)
    for o in obj:
        if o[0]:
            objects.append(o[1].as_ufl())
        else:
            objects.append(o[1])
    return objects
