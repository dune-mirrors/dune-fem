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
# from . import plotting

from warnings import warn
class plotting:
    @staticmethod
    def plotPointData(*args,**kwarg):
        warn("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first", DeprecationWarning, stacklevel=2)
        from dune.fem.plotting import plotPointData
        plotPointData(*args,**kwarg)
    def plotComponents(*args,**kwarg):
        warn("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first", DeprecationWarning, stacklevel=2)
        from dune.fem.plotting import plotComponents
        plotComponents(*args,**kwarg)
    def triangulationOfNetwork(*args,**kwarg):
        warn("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first", DeprecationWarning, stacklevel=2)
        from dune.fem.plotting import triangulationOfNetwork
        triangulationOfNetwork(*args,**kwarg)
    def mayaviPointData(*args,**kwarg):
        warn("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first", DeprecationWarning, stacklevel=2)
        from dune.fem.plotting import mayaviPointData
        mayaviPointData(*args,**kwarg)

# finalization of fem module (i.e. calling PETSc finalize etc)
# atexit.register( _fem.__finalizeFemModule__ )

registry = {}

registry["view"] = {
         "adaptive"   : view.adaptiveLeafGridView,
         "filtered"   : view.filteredGridView,
         "geometry"   : view.geometryGridView
     }
registry["space"] = {
         "lagrange"          : space.lagrange,
         "lagrangehp"        : space.lagrangehp,
         "dgonb"             : space.dgonb,
         "dgonbhp"           : space.dgonbhp,
         "dglegendre"        : space.dglegendre,
         "dglegendrehp"      : space.dglegendrehp,
         "dglagrange"        : space.dglagrange,
         "dglagrangelobatto" : space.dglagrangelobatto,
         "finitevolume"      : space.finiteVolume,
         "p1bubble"          : space.p1Bubble,
         "combined"          : space.combined,
         "product"           : space.product,
         "bdm"               : space.bdm,
         "raviartthomas"     : space.raviartThomas,
         "rannacherTurek"    : space.rannacherTurek
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

def assemble(form,space=None,gridView=None,order=None):
    import ufl
    from ufl.log import UFLException
    from ufl.equation import Equation
    from ufl.algorithms.analysis import extract_arguments_and_coefficients
    from ufl.algorithms.estimate_degrees import estimate_total_polynomial_degree
    try:
        params = [*form[1:]]
        form = form[0]
    except TypeError: # UFLException is thrown if an expression instead of a form is provided
        params = []
        pass
    if isinstance(form, Equation):
        functional = form.rhs
        form = form.lhs
    else:
        functional = None

    args, cc = extract_arguments_and_coefficients(form)
    arity = len(args)
    assert arity == 0 or arity == 1 or arity == 2
    assert functional is None or arity == 2

    if arity == 0:
        # woud be nice to extend this to multiple integrals
        # The whole thing could possibly be implemented using a galerkin
        # operator with test functions in a FV space? This would be similar
        # to the arity 1 case but would then include summing up the result vector.
        # Efficiency is the issue but it might provide a way of including
        # integrals over boundaries or possibly internal boundaries. Or we
        # define a space with a single constant 1 basis function.
        assert len(form.integrals()) == 1, "can only integrate forms with single integral"
        if gridView is None:
            assert len(cc) > 0, "to integrate an expression  a 'gridView' has to be provided or the form must contain a grid function."
            gridView = cc[0].gridView
        if order is None:
            ### assert len(cc) > 0, "to integrate an expression not containing a grid functiona quadrature 'order' has to be provided."
            ### order = max( gf.order for gf in cc if hasattr(gf,"order") )
            # the ufl estimate is quite high,
            # e.g. degree(x)=1, degree(x.x)=2, degree(sin(x))=3=degree(x)+2, degree(sin(x.x))=4=degree(x.x)+2
            # so degree( (df-exact)**2 ) = 8 if exact=sin(x.x) independent of degree(df).
            # With for example df linear Lagrange I would assume that integrating with order=5 is more than enough.
            # If exact=xy and df is linear (exact-df)**2*dx would be integrated with order 4
            # - I would probably have taken 5 again which have been stupid...
            order = estimate_total_polynomial_degree( form )
        return dune.fem.function.integrate(gridView, form.integrals()[0].integrand(), order=order)
    else:
        v = args[0]
        if not space:
            try:
                space = v.ufl_function_space()
            except:
                raise AttributeError("can not access space from given form")
        if arity == 1:
            # todo: implement this on the C++ side - we use a Galerkin operator as a stopgap solution
            u = ufl.TrialFunction(space)
            op = dune.fem.operator.galerkin( [u*v*ufl.dx == form] + params )
            b = space.zero.copy()
            op(space.zero,b)
            return b
        else:
            if functional is not None:
                op = dune.fem.operator.galerkin( [form == functional] + params )
                A = dune.fem.operator.linear(op)
                op.jacobian(space.zero,A)
                b = space.zero.copy()
                op(space.zero,b)
                b *= -1
                return A,b
            else:
                op = dune.fem.operator.galerkin( [form] + params )
                A = dune.fem.operator.linear(op)
                op.jacobian(space.zero,A)
                return A
