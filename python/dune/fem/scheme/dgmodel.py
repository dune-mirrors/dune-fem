from __future__ import division, print_function, unicode_literals

from dune.ufl.codegen import generateMethod
from ufl import SpatialCoordinate, Coefficient, replace, diff, as_vector
from ufl.core.expr import Expr
from dune.source.cplusplus import Variable, UnformattedExpression, AccessModifier
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices, expand_derivatives

def codeDG(self):
    code = self._code()

    u = self.trialFunction
    ubar = Coefficient(u.ufl_function_space())
    penalty = self.penalty
    if penalty is None:
        penalty = 1
    if isinstance(penalty,Expr):
        if penalty.ufl_shape == ():
            penalty = as_vector([penalty])
        try:
            penalty = expand_indices(expand_derivatives(expand_compounds(penalty)))
        except:
            pass
        assert penalty.ufl_shape == u.ufl_shape
        dmPenalty = as_vector([
                    replace(
                      expand_derivatives( diff(replace(penalty,{u:ubar}),ubar) )[i,i],
                    {ubar:u} )
                 for i in range(u.ufl_shape[0]) ])
    else:
        dmPenalty = None

    code.append(AccessModifier("public"))
    x = SpatialCoordinate(self.space.cell())
    predefined = {}
    self.predefineCoefficients(predefined,x)
    spatial = Variable('const auto', 'y')
    predefined.update( {x: UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( x ) )') })
    generateMethod(code, penalty,
            'RRangeType', 'penalty',
            args=['const Point &x',
                  'const DRangeType &u'],
            targs=['class Point','class DRangeType'], static=False, const=True,
            predefined=predefined)
    generateMethod(code, dmPenalty,
            'RRangeType', 'linPenalty',
            args=['const Point &x',
                  'const DRangeType &u'],
            targs=['class Point','class DRangeType'], static=False, const=True,
            predefined=predefined)
    return code

def transform(space,penalty):
    def transform_(model):
        if model.baseName == "modelDG":
            return
        model._code = model.code
        model.code  = lambda *args,**kwargs: codeDG(*args,**kwargs)
        model.space = space
        model.penalty = penalty
        model.baseName = "modelDG"
        model.modelWrapper = "DGDiffusionModelWrapper< Model >"
        # model.baseSignature = []
        # if penalty is not None:
        #     model.baseSignature += [penalty]
    return [transform_,[penalty]]
