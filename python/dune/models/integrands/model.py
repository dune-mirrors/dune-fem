from __future__ import division, print_function, unicode_literals

from dune.ufl import codegen
from dune.source.cplusplus import Include, Method, TypeAlias, Variable, Declaration
from dune.source.cplusplus import assign, construct, coordinate, dereference, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_
from dune.source.fem import fieldTensorType
from dune.source.algorithm.extractincludes import extractIncludesFromStatements


class Integrands(codegen.ModelClass):
    def __init__(self, trialFunction, domainValue, rangeValue=None, uflExpr=None, virtualize=True):
        """construct new integrands

        Args:
            domainValue:      structure of domain value tuple
            rangeVlue:        structure of range value tuple

        Returns:
            Integrands: newly constructed integrands

        The tuples domainValue and rangeValue contain the shapes of the
        corresponding value types for these integrands.

        If parameter names are provided, the tuple must have the same
        length as the constants tuple. For constants without parameter name,
        pass None instead of a string.
        """
        codegen.ModelClass.__init__(self, "Integrands", uflExpr, virtualize)

        self.form = uflExpr[0]

        if rangeValue is None:
            rangeValue = domainValue

        domainValue = tuple(domainValue)
        rangeValue  = tuple(rangeValue)
        self.trialFunction = trialFunction
        self.domainValueTuple = 'std::tuple< ' + ', '.join(fieldTensorType(v) for v in domainValue) + ' >'
        self.rangeValueTuple = 'std::tuple< ' + ', '.join(fieldTensorType(v) for v in rangeValue) + ' >'
        self.dimDomain = domainValue[0][0]
        self.dimRange = rangeValue[0][0]

        self.field = "double"

        self.interior = None
        self.linearizedInterior = None
        self.boundary = None
        self.linearizedBoundary = None
        self.skeleton = None
        self.linearizedSkeleton = None
        self.nonlinear = True
        self.symmetric = False

        # Added for dirichlet treatment (same as conservationlaw model)
        self.hasDirichletBoundary = False
        self.hasNeumanBoundary = False
        self.isDirichletIntersection = None # [return_(False)]
        self.dirichlet = None
        self.dDirichlet = None
        self.arg_i = Variable("const IntersectionType &", "intersection")
        self.arg_bndId = Variable("int", "bndId")
        self.arg_r = Variable("RRangeType &", "result")
        self.arg_dr = Variable("RJacobianRangeType &", "result")
        self.arg_x = Variable("const Point &", "x")

        self.baseName = 'integrands'
        self.baseSignature = []

    def signature(self):
        return self.form.signature()

    def methods(self,code):
        code.append(TypeAlias("DomainValueType", self.domainValueTuple))
        code.append(TypeAlias("RangeValueType", self.rangeValueTuple))
        code.append(Declaration(Variable("bool", "_nonlinear"), initializer=self.nonlinear, static=True, constexpr=True))
        code.append(Method('bool', 'nonlinear', args=[], code=["return _nonlinear;"], const=True))

        if self.interior is not None:
            code.append(Method('RangeValueType', 'interior', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.interior, const=True))
            code.append(Method('auto', 'linearizedInterior', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.linearizedInterior, const=True))

        if self.boundary is not None:
            code.append(Method('RangeValueType', 'boundary', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.boundary, const=True))
            code.append(Method('auto', 'linearizedBoundary', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.linearizedBoundary, const=True))

        if self.skeleton is not None:
            code.append(Method('std::pair< RangeValueType, RangeValueType >', 'skeleton', targs=['class Point'], args=['const Point &xIn', 'const DomainValueType &uIn', 'const Point &xOut', 'const DomainValueType &uOut'], code=self.skeleton, const=True))
            code.append(Method('auto', 'linearizedSkeleton', targs=['class Point'], args=['const Point &xIn', 'const DomainValueType &uIn', 'const Point &xOut', 'const DomainValueType &uOut'], code=self.linearizedSkeleton, const=True))

        # added for dirichlet treatment - same as conservationlaw model
        if self.hasDirichletBoundary is not None:
            code.append(TypeAlias("RRangeType",'Dune::FieldVector< double, '+ str(self.dimRange) + ' > '))
            code.append(TypeAlias("RJacobianRangeType",'Dune::FieldMatrix< double, ' + str(self.dimRange) + ', GridPartType::dimension > '))
            code.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
            # code.append(TypeAlias("BoundaryIdProviderType",\ "Dune::Fem::BoundaryIdGetter< typename GridPartType::GridType >"))
            # idGetter = Variable('BoundaryIdProviderType', "boundaryIdGetter_")
            # code.append(Declaration(idGetter))
            code.append(TypeAlias("DirichletComponentType","std::array<int,"+str(self.dimRange)+">"))
            code.append(Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)))
            code.append(Method('bool', 'isDirichletIntersection', args=[self.arg_i, 'DirichletComponentType &dirichletComponent'],
                code=self.isDirichletIntersection\
                          if self.isDirichletIntersection is not None else\
                     [return_(False)],
                const=True))
            code.append(Method('void', 'dirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_r], code=self.dirichlet, const=True))
            code.append(Method('void', 'dDirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_dr], code=self.dDirichlet, const=True))

    def includes(self):
        incs = set.union(*[extractIncludesFromStatements(stmts) for stmts in (self.interior, self.linearizedInterior, self.boundary, self.linearizedBoundary, self.skeleton, self.linearizedSkeleton)])
        ret = [Include(i) for i in incs]
        ret.sort()
        return ret
