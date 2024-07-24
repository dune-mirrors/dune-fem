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

from dune.fem.deprecated import deprecated

class plotting:
    @staticmethod
    def plotPointData(*args,**kwarg):
        deprecated("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first")
        from dune.fem.plotting import plotPointData
        plotPointData(*args,**kwarg)
    def plotComponents(*args,**kwarg):
        deprecated("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first")
        from dune.fem.plotting import plotComponents
        plotComponents(*args,**kwarg)
    def triangulationOfNetwork(*args,**kwarg):
        deprecated("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first")
        from dune.fem.plotting import triangulationOfNetwork
        triangulationOfNetwork(*args,**kwarg)
    def mayaviPointData(*args,**kwarg):
        deprecated("This use of dune.fem.plotting is deprecated. Import dune.fem.plotting first")
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
         "composite"         : space.composite,
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
         # "global"     : function.globalFunction,
         # "local"      : function.localFunction,
         # "cpp"        : function.cppFunction,
         # "ufl"        : function.uflFunction,
         "gridFunction" : function.gridFunction,
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
    from dune.fem.function._functions import gridFunction
    order = 5 # needs to be derived from f
    try:
        gf = gridFunction(expr=f,gridView=grid, name="tmp", order=order).gf
    except Exception as e: # AttributeError as e:
        gf = None
    return gf
_writeVTKDispatcher.append(vtkDispatchUFL)

def assemble(form,space=None,gridView=None,order=None):
    from dune.ufl import DirichletBC
    import ufl
    from ufl.equation import Equation
    from ufl.algorithms.analysis import extract_arguments_and_coefficients
    from ufl.algorithms.estimate_degrees import estimate_total_polynomial_degree
    if type(form)==tuple or type(form)==list:
        params = [*form[1:]]
        form = form[0]
        if not (isinstance(form,Equation) or isinstance(form,ufl.Form)):
            raise ValueError(
"""If the first argument is a list then the first entry must be either a ufl Form or an Equation"""
            )
        if isinstance(form, ufl.core.expr.Expr):
            raise ValueError("must provide a form or equation not a ufl expression - forgot to multiply with ufl.dx?")
        if not all([isinstance(x,DirichletBC) for x in params]):
            raise ValueError(
"""if the first argument is a list or tuple it must be of the form [ufl_form,DirichletBC,...,DirichletBC]"""
            )
    else:
        params = []
        if isinstance(form, ufl.core.expr.Expr):
            raise ValueError("must provide a form or equation not a ufl expression - forgot to multiply with ufl.dx?")
        if not (isinstance(form,Equation) or isinstance(form,ufl.Form)):
            raise ValueError(
"""If the first argument is a list then the first entry must be either a ufl Form or an Equation"""
            )
        pass

    if isinstance(form,Equation):
        form = form.lhs - form.rhs
        wasEqn = True
    else:
        wasEqn = False
    args, cc = extract_arguments_and_coefficients(form)
    arity = len(args)
    if wasEqn and arity<2:
        raise ValueError("An equation must have arity 2, i.e., contain a bilinear form")

    if arity == 0:
        # woud be nice to extend this to multiple integrals
        # The whole thing could possibly be implemented using a galerkin
        # operator with test functions in a FV space? This would be similar
        # to the arity 1 case but would then include summing up the result vector.
        # Efficiency is the issue but it might provide a way of including
        # integrals over boundaries or possibly internal boundaries. Or we
        # define a space with a single constant 1 basis function.

        if len(form.integrals()) > 1:
            raise ValueError("currently only forms with single integral can be integrated")
        if gridView is None:
            try:
                gridView = cc[0].gridView
            except TypeError:
                raise ValueError(
"""to integrate an expression a 'gridView' has to be provided or the form must contain a grid function."""
                )
        if order is None:
            # the ufl estimate is quite high,
            # e.g. degree(x)=1, degree(x.x)=2, degree(sin(x))=3=degree(x)+2, degree(sin(x.x))=4=degree(x.x)+2
            # so degree( (df-exact)**2 ) = 8 if exact=sin(x.x) independent of degree(df).
            # With for example df linear Lagrange I would assume that integrating with order=5 is more than enough.
            order = estimate_total_polynomial_degree( form )
        return integrate(form.integrals()[0].integrand(), gridView=gridView, order=order)
    else:
        v = args[0]
        if not space:
            try:
                space = v.ufl_function_space()
            except:
                raise ValueError("can not access space from given form - provide a space as argument")
        if arity == 1:
            # todo: implement this on the C++ side - we use a Galerkin operator as a stopgap solution
            u = ufl.TrialFunction(space) # this is not good - the space might not be available
            op = dune.fem.operator.galerkin( [ufl.inner(u,v)*ufl.dx - form] + params )
            b = space.zero.copy()
            op(space.zero,b)
            b *= -1 # note: using u*v*ufl.dx + form to avoid the *=-1 fails since the boundary values would have the wrong sign
            return b
        else:
            op = dune.fem.operator.galerkin( [form] + params )
            A = op.linear()
            zeroU = args[1].ufl_function_space().zero
            if not wasEqn:
                computeRHS = not (form==0).rhs == 0
                # computeRHS = not compute_form_rhs(form).empty()
            else:
                computeRHS = True
            if computeRHS:
                b = space.zero.copy()
                op.jacobian(zeroU,A,b)
                return A,b
            else:
                op.jacobian(zeroU,A)
                return A

def integrate(expr, gridView=None, order=None):
    from dune.fem.function._functions import _integrate
    return _integrate(gridView, expr, order)
