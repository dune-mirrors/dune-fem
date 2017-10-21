from __future__ import division, print_function

from types import MethodType, ModuleType

from ufl import Coefficient, FacetNormal, Form, SpatialCoordinate
from ufl import CellVolume, MinCellEdgeLength, MaxCellEdgeLength
from ufl import FacetArea, MinFacetEdgeLength, MaxFacetEdgeLength
from ufl import action, derivative
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.replace import Replacer
from ufl.constantvalue import IntValue, Zero
from ufl.corealg.map_dag import map_expr_dags
from ufl.equation import Equation
from ufl.differentiation import Grad

from dune.source.builtin import get, hybridForEach, make_pair, make_index_sequence, make_shared
from dune.source.cplusplus import AccessModifier, Declaration, Constructor, EnumClass, Include, InitializerList, Method, NameSpace, Struct, TypeAlias, UnformattedExpression, Using, Variable
from dune.source.cplusplus import assign, construct, coordinate, dereference, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_
from dune.source.cplusplus import SourceWriter
from dune.source.algorithm.extractincludes import extractIncludesFromStatements
from dune.source.algorithm.extractvariables import extractVariablesFromExpressions, extractVariablesFromStatements

from dune.ufl import codegen
from dune.ufl.gatherderivatives import gatherDerivatives
from dune.ufl.linear import splitForm
import dune.ufl.tensors as tensors

class Integrands():
    def __init__(self, signature, domainValue, rangeValue=None):
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
        self._constants = []
        self._coefficients = []
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

    def addCoefficient(self, dimRange, field="double"):
        idx = len(self._coefficients)
        self._coefficients.append({'dimRange': dimRange, 'field': field})
        return idx

    def addConstant(self, cppType):
        idx = len(self._constants)
        self._constants.append(cppType)
        return idx

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
            coefficientSpaces = [('Dune::Fem::FunctionSpace< double,' + SourceWriter.cpp_fields(c['field']) + ', GlobalCoordinateType::dimension, ' + str(c['dimRange']) + ' >') for c in self._coefficients]
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


def generateCode(predefined, testFunctions, tensorMap, tempVars=True):
    # build list of all expressions to compile
    expressions = []
    for phi in testFunctions:
        if (phi,) not in tensorMap:
            continue
        tensor = tensorMap[(phi,)]
        keys = tensor.keys()
        expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)

    # extract generated code for expressions and build values
    values = []
    for phi in testFunctions:
        value = tensors.fill(phi.ufl_shape, makeExpression(0))
        if (phi,) in tensorMap:
            tensor = tensorMap[(phi,)]
            keys = tensor.keys()
            for i, r in zip(keys, results[:len(keys)]):
                value = tensors.setItem(value, i, r)
            results = results[len(keys):]
        values += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values


def generateLinearizedCode(predefined, testFunctions, trialFunctionMap, tensorMap, tempVars=True):
    """generate code for a bilinear form

    Args:
        predefined:       list of predefined arguments or coefficients
        testFunctions:    list of arguments to interpret as test functions
        trialFunctionMap: map of variable to list of arguments to interpret as trial functions
        tensorMap:        map of expression tensors of shape (testFunction x trialFunction)
        tempVars:         introduce temporary variables during code generation
    """

    # build list of all expressions to compile
    expressions = []
    for var, trialFunctions in trialFunctionMap.items():
        for phi in testFunctions:
            for psi in trialFunctions:
                if (phi, psi) not in tensorMap:
                    continue
                tensor = tensorMap[(phi, psi)]
                keys = tensor.keys()
                expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, tempVars=tempVars)

    # extract generated code for expressions and build values
    values = {}
    for var, trialFunctions in trialFunctionMap.items():
        values[var] = []
        for phi in testFunctions:
            value = tensors.fill(phi.ufl_shape, None)
            for idx in range(len(trialFunctions)):
                psi = trialFunctions[idx]
                if (phi, psi) in tensorMap:
                    tensor = tensorMap[(phi, psi)]
                    keys = tensor.keys()
                    for ij, r in zip(keys, results[:len(keys)]):
                        if isinstance(tensor[ij], Zero):
                            continue
                        i = ij[:len(phi.ufl_shape)]
                        j = ij[len(phi.ufl_shape):]
                        if isinstance(tensor[ij], IntValue) and int(tensor[ij]) == 1:
                            r = var[idx][j]
                        else:
                            r = r * var[idx][j]
                        s = tensors.getItem(value, i)
                        s = r if s is None else s + r
                        value = tensors.setItem(value, i, s)
                    results = results[len(keys):]
            value = tensors.apply(lambda v : makeExpression(0) if v is None else v, phi.ufl_shape, value)
            values[var] += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values


def generateUnaryCode(predefined, testFunctions, tensorMap, tempVars=True):
    preamble, values = generateCode(predefined, testFunctions, tensorMap, tempVars=tempVars)
    return preamble + [return_(construct('RangeValueType', *values, brace=True))]


def generateUnaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, tempVars=True):
    if tensorMap is None:
        return [return_(lambda_(args=['const DomainValueType &phi'], code=return_(construct('RangeValueType', *[0 for i in range(len(testFunctions))], brace=True))))]

    var = Variable('std::tuple< RangeType, JacobianRangeType >', 'phi')
    preamble, values = generateLinearizedCode(predefined, testFunctions, {var: trialFunctions}, tensorMap, tempVars=tempVars)
    capture = extractVariablesFromExpressions(values[var]) - {var}
    return preamble + [return_(lambda_(capture=capture, args=['const DomainValueType &phi'], code=return_(construct('RangeValueType', *values[var], brace=True))))]


def generateBinaryCode(predefined, testFunctions, tensorMap, tempVars=True):
    restrictedTestFunctions = [phi('+') for phi in testFunctions] + [phi('-') for phi in testFunctions]
    preamble, values = generateCode(predefined, restrictedTestFunctions, tensorMap, tempVars=tempVars)
    return preamble + [return_(make_pair(construct('RangeValueType', *values[:len(testFunctions)], brace=True), construct('RangeValueType', *values[len(testFunctions):], brace=True)))]


def generateBinaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, tempVars=True):
    restrictedTestFunctions = [phi('+') for phi in testFunctions] + [phi('-') for phi in testFunctions]

    trialFunctionsIn = [psi('+') for psi in trialFunctions]
    trialFunctionsOut = [psi('-') for psi in trialFunctions]

    if tensorMap is None:
        value = construct('RangeValueType', *[0 for i in range(len(testFunctions))], brace=True)
        tensorIn = lambda_(args=['const DomainValueType &phiIn'], code=return_(make_pair(value, value)))
        tensorOut = lambda_(args=['const DomainValueType &phiOut'], code=return_(make_pair(value, value)))
        return [return_(make_pair(tensorIn, tensorOut))]

    varIn = Variable('std::tuple< RangeType, JacobianRangeType >', 'phiIn')
    varOut = Variable('std::tuple< RangeType, JacobianRangeType >', 'phiOut')
    preamble, values = generateLinearizedCode(predefined, restrictedTestFunctions, {varIn: trialFunctionsIn, varOut: trialFunctionsOut}, tensorMap, tempVars=tempVars)

    captureIn = extractVariablesFromExpressions(values[varIn]) - {varIn}
    captureOut = extractVariablesFromExpressions(values[varOut]) - {varOut}

    tensorIn = lambda_(capture=captureIn, args=['const DomainValueType &phiIn'], code=return_(make_pair(construct('RangeValueType', *values[varIn][:len(testFunctions)], brace=True), construct('RangeValueType', *values[varIn][len(testFunctions):], brace=True))))
    tensorOut = lambda_(capture=captureOut, args=['const DomainValueType &phiOut'], code=return_(make_pair(construct('RangeValueType', *values[varOut][:len(testFunctions)], brace=True), construct('RangeValueType', *values[varOut][len(testFunctions):], brace=True))))

    return preamble + [return_(make_pair(tensorIn, tensorOut))]


def compileUFL(equation, tempVars=True):
    form = equation.lhs - equation.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    x = SpatialCoordinate(form.ufl_cell())
    n = FacetNormal(form.ufl_cell())

    cellVolume = CellVolume(form.ufl_cell())
    maxCellEdgeLength = MaxCellEdgeLength(form.ufl_cell())
    minCellEdgeLength = MinCellEdgeLength(form.ufl_cell())

    facetArea = FacetArea(form.ufl_cell())
    maxFacetEdgeLength = MaxFacetEdgeLength(form.ufl_cell())
    minFacetEdgeLength = MinFacetEdgeLength(form.ufl_cell())

    phi = form.arguments()[0]

    u = form.arguments()[1]
    ubar = Coefficient(u.ufl_function_space())

    derivatives = gatherDerivatives(form, [phi, u])

    derivatives_phi = derivatives[0]
    derivatives_u = derivatives[1]
    derivatives_ubar = map_expr_dags(Replacer({u: ubar}), derivatives_u)

    integrands = Integrands(form.signature(), (d.ufl_shape for d in derivatives_u), (d.ufl_shape for d in derivatives_phi))
    try:
        integrands.field = u.ufl_function_space().field()
    except AttributeError:
        pass

    constants = dict()
    coefficients = dict()
    for coefficient in set(form.coefficients()):
        if coefficient.is_cellwise_constant():
            dimRange = (1 if len(coefficient.ufl_shape) == 0 else coefficient.ufl_shape[0])
            constants[coefficient] = integrands.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >')
        else:
            try:
                coefficients[coefficient] = integrands.addCoefficient(coefficient.ufl_shape[0], coefficient.ufl_function_space().field())
            except AttributeError:
                coefficients[coefficient] = integrands.addCoefficient(coefficient.ufl_shape[0])

    integrals = splitForm(form, [phi])

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))
    linearizedIntegrals = splitForm(dform, [phi, u])

    if not set(integrals.keys()) <= {'cell', 'exterior_facet', 'interior_facet'}:
        raise Exception('unknown integral encountered in ' + str(set(integrals.keys())) + '.')

    def predefineCoefficients(predefined, x, side=None):
        for coefficient, idx in coefficients.items():
            for derivative in integrands.coefficient(idx, x, side=side):
                if side is None:
                    predefined[coefficient] = derivative
                elif side == 'Side::in':
                    predefined[coefficient('+')] = derivative
                elif side == 'Side::out':
                    predefined[coefficient('-')] = derivative
                coefficient = Grad(coefficient)

    if 'cell' in integrals.keys():
        arg = Variable(integrands.domainValueTuple(), 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.interior = generateUnaryCode(predefined, derivatives_phi, integrals['cell'], tempVars=tempVars)

        predefined = {derivatives_ubar[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.linearizedInterior = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('cell'), tempVars=tempVars)

    if 'exterior_facet' in integrals.keys():
        arg = Variable(integrands.domainValueTuple(), 'u')

        predefined = {derivatives_u[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[n] = integrands.facetNormal('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.boundary = generateUnaryCode(predefined, derivatives_phi, integrals['exterior_facet'], tempVars=tempVars);

        predefined = {derivatives_ubar[i]: arg[i] for i in range(len(derivatives_u))}
        predefined[x] = integrands.spatialCoordinate('x')
        predefined[n] = integrands.facetNormal('x')
        predefined[cellVolume] = integrands.cellVolume()
        predefined[maxCellEdgeLength] = maxEdgeLength(integrands.cellGeometry())
        predefined[minCellEdgeLength] = minEdgeLength(integrands.cellGeometry())
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'x')
        integrands.linearizedBoundary = generateUnaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('exterior_facet'), tempVars=tempVars)

    if 'interior_facet' in integrals.keys():
        argIn = Variable(integrands.domainValueTuple(), 'uIn')
        argOut = Variable(integrands.domainValueTuple(), 'uOut')

        predefined = {derivatives_u[i](s): arg[i] for i in range(len(derivatives_u)) for s, arg in (('+', argIn), ('-', argOut))}
        predefined[x] = integrands.spatialCoordinate('xIn')
        predefined[n('+')] = integrands.facetNormal('xIn')
        predefined[cellVolume('+')] = integrands.cellVolume('Side::in')
        predefined[cellVolume('-')] = integrands.cellVolume('Side::out')
        predefined[maxCellEdgeLength('+')] = maxEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[maxCellEdgeLength('-')] = maxEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[minCellEdgeLength('+')] = minEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[minCellEdgeLength('-')] = minEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'xIn', 'Side::in')
        predefineCoefficients(predefined, 'xOut', 'Side::out')
        integrands.skeleton = generateBinaryCode(predefined, derivatives_phi, integrals['interior_facet'], tempVars=tempVars)

        predefined = {derivatives_ubar[i](s): arg[i] for i in range(len(derivatives_u)) for s, arg in (('+', argIn), ('-', argOut))}
        predefined[x] = integrands.spatialCoordinate('xIn')
        predefined[n('+')] = integrands.facetNormal('xIn')
        predefined[cellVolume('+')] = integrands.cellVolume('Side::in')
        predefined[cellVolume('-')] = integrands.cellVolume('Side::out')
        predefined[maxCellEdgeLength('+')] = maxEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[maxCellEdgeLength('-')] = maxEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[minCellEdgeLength('+')] = minEdgeLength(integrands.cellGeometry('Side::in'))
        predefined[minCellEdgeLength('-')] = minEdgeLength(integrands.cellGeometry('Side::out'))
        predefined[facetArea] = integrands.facetArea()
        predefined[maxFacetEdgeLength] = maxEdgeLength(integrands.facetGeometry())
        predefined[minFacetEdgeLength] = minEdgeLength(integrands.facetGeometry())
        predefined.update({c: integrands.constant(i) for c, i in constants.items()})
        predefineCoefficients(predefined, 'xIn', 'Side::in')
        predefineCoefficients(predefined, 'xOut', 'Side::out')
        integrands.linearizedSkeleton = generateBinaryLinearizedCode(predefined, derivatives_phi, derivatives_u, linearizedIntegrals.get('interior_facet'), tempVars=tempVars)

    coefficients.update(constants)
    return integrands, coefficients


def setConstant(integrands, index, value):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setConstant(index, value)


def setCoefficient(integrands, index, coefficient):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setCoefficient(index, coefficient)


def load(grid, integrands, renumbering=None, tempVars=True):
    from dune.common.hashit import hashIt

    if isinstance(integrands, Equation):
        integrands, renumbering = compileUFL(integrands, tempVars=tempVars)

    name = 'integrands_' + integrands.signature + '_' + hashIt(grid._typeName)

    code = [Include('config.h')]
    code += [Include(i) for i in grid._includes]
    #code.append(Include("dune/fem/misc/boundaryidprovider.hh"))

    code += integrands.includes()
    code.append(Include("dune/python/pybind11/pybind11.h"))
    code.append(Include("dune/python/pybind11/extensions.h"))
    code.append(Include("dune/fempy/py/grid/gridpart.hh"))

    if integrands._coefficients:
        code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
    code.append(Include("dune/fempy/py/integrands.hh"))

    nameSpace = NameSpace('Integrands_' + integrands.signature)
    nameSpace.append(integrands.code())
    code.append(nameSpace)

    code.append(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + grid._typeName + ' >'))
    if integrands._coefficients:
        rangeTypes = ['Dune::FieldVector< ' + SourceWriter.cpp_fields(c['field']) + ', ' + str(c['dimRange']) + ' >' for c in integrands._coefficients]
        coefficients = ['Dune::FemPy::VirtualizedGridFunction< GridPart, ' + r + ' >' for r in rangeTypes]
    else:
        coefficients = []
    code.append(TypeAlias('Integrands', nameSpace.name + '::Integrands< ' + ', '.join(['GridPart'] + coefficients) + ' >'))

    writer = SourceWriter()
    writer.emit(code);

    writer.openPythonModule(name)
    writer.emit('auto cls = Dune::FemPy::registerIntegrands< Integrands >( module );')
    # TODO
    writer.emit('cls.def( "__init__", [] ( Integrands &self ) { new (&self) Integrands(); } );')
    writer.closePythonModule(name)

    source = writer.writer.getvalue()
    writer.close()

    from dune.generator import builder
    module = builder.load(name, source, "integrands")
    setattr(module.Integrands, "_domainValueType", integrands.domainValueTuple())
    setattr(module.Integrands, "_rangeValueType", integrands.rangeValueTuple())
    if (renumbering is not None) and not hasattr(module.Integrands, "_renumbering"):
        module.Integrands._setConstant = module.Integrands.__dict__['setConstant']
        module.Integrands._setCoefficient = module.Integrands.__dict__['setCoefficient']
        #setattr(module.Integrands, "_setConstant", getattr(module.Integrands, "setConstant"))
        #setattr(module.Integrands, "_setCoefficient", getattr(module.Integrands, "setCoefficient"))
        setattr(module.Integrands, '_renumbering', renumbering)
        setattr(module.Integrands, 'setConstant', setConstant)
        setattr(module.Integrands, 'setCoefficient', setCoefficient)
    return module


def create(grid, integrands, renumbering=None, tempVars=True):
    return load(grid, integrands, renumbering=renumbering, tempVars=tempVars).Integrands()
