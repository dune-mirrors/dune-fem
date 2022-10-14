from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import logging
logger = logging.getLogger(__name__)

import dune.common.checkconfiguration as checkconfiguration
from dune.ufl import DirichletBC

def warnNewCartesianIds(view,*args):
    try:
        needsChecking = view.hierarchicalGrid._cartesianConstructionWithIds
    except AttributeError:
        needsChecking = False

    if needsChecking:
        dirichletBCs = [arg for arg in args
                        if isinstance(arg, DirichletBC)]
        if ( len(dirichletBCs) == 1
                 and isinstance(dirichletBCs[0].subDomain,int)
                 and dirichletBCs[0].subDomain==1 ):
            print(
"""
Warning: note that the boundary ids of the Cartesian domain have changed
for some unstructured grids (alu and mmesh). They are now set according to
the convention for structured grids based on the dune reference element, i.e.,
in 2d left=1, right=2, bottom=3, top=4. How to get rid of this warning:

- You can switch to the old ids by passing 'boundary=False' to the
  'cartesianDomain' function. But this option will be removed soon.
- If you want to set boundary conditions on the whole boundary you can remove the
  final parameter ('subDomain') from the 'DirichletBC' class.
- If you are already using the new boundary ids and want only Dirichlet
  conditions on the left side (so the warning is wrongly provided) then
  please add a second 'DirichletBC' instance to the scheme constructor with
  a boundary id above 2^dim - sorry about this...
""")

def conservationlaw(view, equation, *args, **kwargs):
    import ufl
    import dune.ufl
    import dune.models.conservationlaw as conservationlaw

    coefficients = kwargs.pop('coefficients', dict())

    warnNewCartesianIds(view,*args)
    Model = conservationlaw.load(view, equation, *args, **kwargs).Model
    # the following needs to be done for the linearized problem:
    # if isinstance(equation, ufl.equation.Equation):
    #     lhs = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(equation.lhs)))
    #     if lhs == ufl.adjoint(lhs):
    #         setattr(Model, 'symmetric', 'true')
    #    else:
    #        setattr(Model, 'symmetric', 'false')

    return Model(coefficients=coefficients)

# deprecated, use conservationlaw instead
def elliptic(view, equation, *args, **kwargs):
    return conservationlaw( view, equation, *args, **kwargs)

def integrands(view, form, *args, **kwargs):
    import ufl
    import dune.ufl
    import dune.models.integrands as integrands

    tempVars   = kwargs.pop('tempVars', True)
    virtualize = kwargs.pop('virtualize',True)
    modelPatch = kwargs.pop('modelPatch',None)
    includes   = kwargs.pop('includes',None)

    warnNewCartesianIds(view,*args)
    Integrands = integrands.load(view, form, *args,
                     tempVars=tempVars,virtualize=virtualize,
                     modelPatch=modelPatch,includes=includes)

    return Integrands(*args, **kwargs)
