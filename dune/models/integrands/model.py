from __future__ import division, print_function, unicode_literals

from dune.source.builtin import get, hybridForEach, make_pair, make_index_sequence, make_shared
from dune.source.cplusplus import AccessModifier, Declaration, Constructor, EnumClass, Include, InitializerList, Method, Struct, TypeAlias, UnformattedExpression, Variable
from dune.source.cplusplus import assign, construct, coordinate, dereference, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_
from dune.source.cplusplus import SourceWriter
from dune.source.algorithm.extractincludes import extractIncludesFromStatements

class Integrands():
    def __init__(self, signature, domainValue, rangeValue=None, constants=None, coefficients=None):
        """construct new integrands

        Args:
            signature:    unique signature for these integrands
            domainValue:  structure of domain value tuple
            rangeVlue:    structure of range value tuple

        Returns:
            Integrands: newly constructed integrands

        The tuples domainValue and rangeValue contain the shapes of the
        corresponding value types for these integrands.
        """
        if rangeValue is None:
            rangeValue = domainValue

        self.signature = signature
        self.domainValue = tuple(domainValue)
        self.rangeValue = tuple(rangeValue)

        self.field = "double"
        self._constants = [] if constants is None else list(constants)
        self._coefficients = [] if coefficients is None else list(coefficients)
        self.init = None
        self.vars = None

        self.interior = None
        self.linearizedInterior = None
        self.boundary = None
        self.linearizedBoundary = None
        self.skeleton = None
        self.linearizedSkeleton = None

        self._derivatives = [('RangeType', 'evaluate'), ('JacobianRangeType', 'jacobian'), ('HessianRangeType', 'hessian')]

    def _cppTensor(self, shape):
        if len(shape) == 1:
            return 'Dune::FieldVector< double, ' + str(shape[0]) + ' >'
        elif len(shape) == 2:
            return 'Dune::FieldMatrix< double, ' + str(shape[0]) + ', ' + str(shape[1]) + ' >'
        elif len(shape) == 3:
            return 'Dune::FieldVector< Dune::FieldMatrix< double, ' + str(shape[1]) + ', ' + str(shape[2]) + ' >, ' + str(shape[0]) + ' >'
        else:
            raise ValueError('No C++ type defined for tensors of shape ' + str(shape) + '.')

    def domainValueTuple(self):
        return 'std::tuple< ' + ', '.join([self._cppTensor(v) for v in self.domainValue]) + ' >'

    def rangeValueTuple(self):
        return 'std::tuple< ' + ', '.join([self._cppTensor(v) for v in self.rangeValue]) + ' >'

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def coefficient(self, idx, x, side=None):
        targs = [str(idx)]
        if side is not None:
            targs.append(side)
        return (UnformattedExpression('typename CoefficientFunctionSpaceType< ' + str(idx) + ' >::' + t, n + 'Coefficient< ' + ', '.join(targs) + ' >( ' + x + ' )') for t, n in self._derivatives)

    def spatialCoordinate(self, x):
        return UnformattedExpression('GlobalCoordinateType', 'entity().geometry().global( Dune::Fem::coordinate( ' + x + ' ) )')

    def facetNormal(self, x):
        return UnformattedExpression('GlobalCoordinateType', 'intersection_.unitOuterNormal( ' + x + '.localPosition() )')

    def cellVolume(self, side=None):
        entity = 'entity()' if side is None else 'entity_[ static_cast< std::size_t >( ' + side + ' ) ]'
        return UnformattedExpression('auto', entity + '.geometry().volume()')

    def cellGeometry(self, side=None):
        entity = 'entity()' if side is None else 'entity_[ static_cast< std::size_t >( ' + side + ' ) ]'
        return UnformattedExpression('auto', entity + '.geometry()')

    def facetArea(self):
        return UnformattedExpression('auto', 'intersection_.geometry().volume()')

    def facetGeometry(self):
        return UnformattedExpression('auto', 'intersection_.geometry()')

    def code(self, name='Integrands', targs=[]):
        code = Struct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']))

        code.append(TypeAlias("GridPartType", "GridPart"))

        code.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))

        code.append(TypeAlias("GlobalCoordinateType", "typename EntityType::Geometry::GlobalCoordinate"))

        code.append(TypeAlias("DomainValueType", self.domainValueTuple()))
        code.append(TypeAlias("RangeValueType", self.rangeValueTuple()))

        constants = ["std::shared_ptr< " + c + " >" for c in self._constants]
        if constants:
            code.append(TypeAlias("ConstantTupleType", "std::tuple< " + ", ".join(constants) + " >"))
            code.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantTupleType >::element_type", targs=["std::size_t i"]))
        else:
            code.append(TypeAlias("ConstantTupleType", "std::tuple<>"))

        if self._coefficients:
            coefficientSpaces = [('Dune::Fem::GridFunctionSpace< GridPartType, ' + c + ' >') for c in self._coefficients]
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " +", ".join(coefficientSpaces) + " >"))
            code.append(TypeAlias("CoefficientTupleType", "std::tuple< Coefficients... >"))

            code.append(TypeAlias("CoefficientFunctionSpaceType", "typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type", targs=["std::size_t i"]))
            for s in ["RangeType", "JacobianRangeType"]:
                code.append(TypeAlias("Coefficient" + s, "typename CoefficientFunctionSpaceType< i >::" + s, targs=["std::size_t i"]))
        else:
            code.append(TypeAlias("CoefficientTupleType", "std::tuple<>"))

        code.append(TypeAlias('CoefficientType', 'typename std::tuple_element< i, CoefficientTupleType >::type', targs=['std::size_t i']))
        code.append(TypeAlias('ConstantType', 'typename std::tuple_element< i, ConstantTupleType >::type::element_type', targs=['std::size_t i']))

        if self.skeleton is not None:
            code.append(EnumClass('Side', ['in = 0u', 'out = 1u'], 'std::size_t'))
            inside = '[ static_cast< std::size_t >( Side::in ) ]'
        else:
            inside = ''

        if self.skeleton is None:
            entity_ = Variable('EntityType', 'entity_')
            insideEntity = entity_
        else:
            entity_ = Variable('std::array< EntityType, 2 >', 'entity_')
            insideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::in )')]
            outsideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::out )')]
        intersection_ = Variable('IntersectionType', 'intersection_')

        constants_ = Variable("ConstantTupleType", "constants_")
        if self.skeleton is None:
            coefficients_ = Variable('std::tuple< typename Coefficients::LocalFunctionType... >', 'coefficients_')
        else:
            coefficients_ = Variable('std::array< std::tuple< typename Coefficients::LocalFunctionType... >, 2 >', 'coefficients_')

        code.append(Constructor(code=hybridForEach(make_index_sequence("std::tuple_size< ConstantTupleType >::value")(), lambda_(args=["auto i"], capture=[Variable('auto', 'this')], code=assign(get("i")(constants_), make_shared("ConstantType< i >")())))))

        entity = Variable('const EntityType &', 'entity')
        initEntity = Method('bool', 'init', args=[entity])
        initEntity.append(assign(insideEntity, entity))
        if self._coefficients:
            initEntity.append('initCoefficients( entity_' + inside + ', coefficients_' + inside + ', std::index_sequence_for< Coefficients... >() );')
        initEntity.append(self.init)
        initEntity.append(return_(True))
        code.append(initEntity)

        intersection = Variable('const IntersectionType &', 'intersection')
        initIntersection = Method('bool', 'init', args=[intersection])
        initIntersection.append(assign(intersection_, intersection))
        if self.skeleton is None:
            initIntersection.append(return_('(intersection.boundary() && init( intersection.inside() ))'))
        else:
            initIntersection.append(assign(insideEntity, UnformattedExpression('EntityType', 'intersection.inside()')))
            if self._coefficients:
                initIntersection.append('initCoefficients( entity_[ static_cast< std::size_t >( Side::in ) ], coefficients_[ static_cast< std::size_t >( Side::in ) ], std::index_sequence_for< Coefficients... >() );')
            initIntersection.append('if( intersection.neighbor() )')
            initIntersection.append('{')
            initIntersection.append('  entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();')
            if self._coefficients:
                initIntersection.append('  initCoefficients( entity_[ static_cast< std::size_t >( Side::out ) ], coefficients_[ static_cast< std::size_t >( Side::out ) ], std::index_sequence_for< Coefficients... >() );')
            initIntersection.append('}')
            initIntersection.append(return_(True))
        code.append(initIntersection)

        if self.interior is not None:
            code.append(Method('RangeValueType', 'interior', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.interior, const=True))
            code.append(Method('auto', 'linearizedInterior', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.linearizedInterior, const=True))

        if self.boundary is not None:
            code.append(Method('RangeValueType', 'boundary', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.boundary, const=True))
            code.append(Method('auto', 'linearizedBoundary', targs=['class Point'], args=['const Point &x', 'const DomainValueType &u'], code=self.linearizedBoundary, const=True))

        if self.skeleton is not None:
            code.append(Method('std::pair< RangeValueType, RangeValueType >', 'skeleton', targs=['class Point'], args=['const Point &xIn', 'const DomainValueType &uIn', 'const Point &xOut', 'const DomainValueType &uOut'], code=self.skeleton, const=True))
            code.append(Method('auto', 'linearizedSkeleton', targs=['class Point'], args=['const Point &xIn', 'const DomainValueType &uIn', 'const Point &xOut', 'const DomainValueType &uOut'], code=self.linearizedSkeleton, const=True))

        code.append(Method('const ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_))), const=True))
        code.append(Method('ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_)))))

        if self._coefficients:
            setCoefficient = Method('void', 'setCoefficient', targs=['std::size_t i'], args=['const CoefficientType< i > &coefficient'])
            if self.skeleton is None:
                setCoefficient.append(assign(get('i')(coefficients_), UnformattedExpression('auto', 'coefficient.localFunction()')))
            else:
                setCoefficient.append(assign(get('i')(coefficients_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::in )')]), UnformattedExpression('auto', 'coefficient.localFunction()')))
                setCoefficient.append(assign(get('i')(coefficients_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::out )')]), UnformattedExpression('auto', 'coefficient.localFunction()')))
            code.append(setCoefficient)

        code.append(Method('const EntityType &', 'entity', const=True, code=return_(insideEntity)))

        code.append(AccessModifier('private'))

        if self._coefficients:
            initCoefficients = Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['const EntityType &entity', 'std::tuple< typename Coefficients::LocalFunctionType... > &coefficients', 'std::index_sequence< i... >'], static=True)
            initCoefficients.append('std::ignore = std::make_tuple( (std::get< i >( coefficients ).init( entity ), i)... );')
            code.append(initCoefficients)

        if self._coefficients:
            for cppType, name in self._derivatives:
                var = Variable('typename CoefficientFunctionSpaceType< i >::' + cppType, 'result')
                if self.skeleton is None:
                    method = Method(var.cppType, name + 'Coefficient', targs=['std::size_t i', 'class Point'], args=['const Point &x'], const=True)
                    method.append(Declaration(var))
                    method.append(UnformattedExpression('void', 'std::get< i >( coefficients_ ).' + name + '( x, ' + var.name + ' );'))
                    method.append(return_(var))
                    code.append(method)
                else:
                    method = Method(var.cppType, name + 'Coefficient', targs=['std::size_t i', 'Side side', 'class Point'], args=['const Point &x'], const=True)
                    method.append(Declaration(var))
                    method.append(UnformattedExpression('void', 'std::get< i >( coefficients_[ static_cast< std::size_t >( side ) ] ).' + name + '( x, ' + var.name + ' )'))
                    method.append(return_(var))
                    code.append(method)

                    method = Method(var.cppType, name + 'Coefficient', targs=['std::size_t i', 'class Point'], args=['const Point &x'], const=True)
                    method.append(return_(UnformattedExpression(var.cppType, name + 'Coefficient< i, Side::in >( x )')))
                    code.append(method)

        code.append(Declaration(entity_), Declaration(intersection_))
        code.append(Declaration(constants_), Declaration(coefficients_))
        if self.vars is not None:
            code += self.vars

        return code

    def includes(self):
        incs = set.union(*[extractIncludesFromStatements(stmts) for stmts in (self.interior, self.linearizedInterior, self.boundary, self.linearizedBoundary, self.skeleton, self.linearizedSkeleton)])
        return [Include(i) for i in incs]
