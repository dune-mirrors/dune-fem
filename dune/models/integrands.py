from __future__ import division, print_function

from ufl import Coefficient, Form
from ufl import action, derivative
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.constantvalue import IntValue, Zero
from ufl.differentiation import Grad

from dune.source.cplusplus import AccessModifier, Declaration, Constructor, EnumClass, InitializerList, Method, Struct, TypeAlias, Variable
from dune.source.cplusplus import construct, lambda_, make_pair, makeExpression, return_
from dune.source.cplusplus import SourceWriter
from dune.source.algorithm.extractvariables import extractVariables

from dune.ufl import codegen
from dune.ufl.linear import splitForm
import dune.ufl.tensors as tensors

class Integrands():
    def __init__(self, dimRange, signature):
        self.dimRange = dimRange
        self.signature = signature
        self.field = "double"
        self.coefficients = []
        self.init = None
        self.vars = None

        self.interior = None
        self.linearizedInterior = None
        self.boundary = None
        self.linearizedBoundary = None
        self.skeleton = None
        self.linearizedSkeleton = None

    def pre(self):
        result = []

        result.append(TypeAlias("GridPartType", "GridPart"))
        result.append(TypeAlias("RangeFieldType", SourceWriter.cpp_fields(self.field)))

        result.append(Declaration(Variable("const int", "dimRange"), str(self.dimRange), static=True))
        result.append(Declaration(Variable("const int", "dimDomain"), "GridPartType::dimensionworld", static=True))
        result.append(Declaration(Variable("const int", "dimLocal"), "GridPartType::dimension", static=True))

        result.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        result.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))
        result.append(TypeAlias("FunctionSpaceType", "Dune::Fem::FunctionSpace< double, RangeFieldType, dimDomain, dimRange >"))

        for s in ["DomainType", "RangeType", "JacobianRangeType"]:
            result.append(TypeAlias(s, "typename FunctionSpaceType::" + s))
        result.append(TypeAlias("ValueType", "std::tuple< RangeType, JacobianRangeType >"))

        if self.coefficients:
            constants = [("std::shared_ptr< Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " > >") for c in self.coefficients if c['constant']]
            result.append(TypeAlias("ConstantsTupleType", "std::tuple< " + ", ".join(constants) + " >"))
            coefficients = [('Dune::Fem::FunctionSpace< double,' + SourceWriter.cpp_fields(c['field']) + ', dimDomain, ' + str(c['dimRange']) + ' >') for c in self.coefficients if not c['constant']]
            result.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " +", ".join(coefficients) + " >"))

            result.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantsTupleType >::element_type", targs=["std::size_t i"]))
            result.append(Declaration(Variable("const std::size_t", "numCoefficients"), "std::tuple_size< CoefficientFunctionSpaceTupleType >::value", static=True))
            result.append(TypeAlias("CoefficientFunctionSpaceType", "typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type", targs=["std::size_t i"]))
            for s in ["RangeType", "JacobianRangeType"]:
                result.append(TypeAlias("Coefficient" + s, "typename CoefficientFunctionSpaceType< i >::" + s, targs=["std::size_t i"]))
        else:
            result.append(Declaration(Variable("const std::size_t", "numCoefficients"), "0u", static=True))
            result.append(TypeAlias("ConstantsTupleType", "std::tuple<>"))

        result.append(TypeAlias('CoefficientType', 'typename std::tuple_element< i, std::tuple< Coefficients... > >::type', targs=['std::size_t i']))
        result.append(TypeAlias('ConstantsType', 'typename std::tuple_element< i, ConstantsTupleType >::type::element_type', targs=['std::size_t i']))

        if self.skeleton is not None:
            result.append(EnumClass('Side', ['in = 0u', 'out = 1u'], 'std::size_t'))
            inside = '[ Side::in ]'
        else:
            inside = ''

        result.append(Constructor(code="constructConstants( std::make_index_sequence< std::tuple_size< ConstantsTupleType >::value >() );"))

        initEntity = Method('bool', 'init', args=['const EntityType &entity'])
        initEntity.append('entity_' + inside + ' = entity;')
        initEntity.append('initCoefficients( entity_' + inside + ', coefficients_' + inside + ', std::make_index_sequence< numCoefficients >() );')
        initEntity.append(self.init)
        initEntity.append(return_(True))
        result.append(initEntity)

        initIntersection = Method('bool', 'init', args=['const IntersectionType &intersection'])
        if self.skeleton is None:
            initIntersection.append(return_('(intersection.boundary() && init( intersection.inside() ))'))
        else:
            initIntersection.append('entity_[ Side::in ] = intersection.inside();')
            initIntersection.append('initCoefficients( entity_[ Side::in ], coefficients_[ Side::in ], std::make_index_sequence< numCoefficients >() );')
            initIntersection.append('if( intersection.neighbor() )')
            initIntersection.append('{')
            initIntersection.append('  entity_[ Side::out ] = intersection.outside();')
            initIntersection.append('  initCoefficients( entity_[ Side::out ], coefficients_[ Side::out ], std::make_index_sequence< numCoefficients >() );')
            initIntersection.append('}')
            initIntersection.append(return_(True))
        result.append(initIntersection)

        return result

    def main(self):
        result = []

        if self.interior is not None:
            result.append(Method('ValueType', 'interior', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.interior, const=True))
            result.append(Method('auto', 'linearizedInterior', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.linearizedInterior, const=True))

        if self.boundary is not None:
            result.append(Method('ValueType', 'boundary', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.boundary, const=True))
            result.append(Method('auto', 'linearizedBoundary', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.linearizedBoundary, const=True))

        if self.skeleton is not None:
            result.append(Method('std::pair< ValueType, ValueType >', 'skeleton', targs=['class Point'], args=['const Point &xIn', 'const ValueType &uIn', 'const Point &xOut', 'const ValueType &uOut'], code=self.skeleton, const=True))
            result.append(Method('auto', 'linearizedSkeleton', targs=['class Point'], args=['const Point &xIn', 'const ValueType &uIn', 'const Point &xOut', 'const ValueType &uOut'], code=self.linearizedSkeleton, const=True))

        return result

    def post(self):
        result = []

        constant = Method('ConstantsType< i > &', 'constant', targs=['std::size_t i'], code=return_('*std::get< i >( constants_ )'))
        result += [constant.variant('const ConstantsType< i > &constant', const=True), constant]

        if self.skeleton is None:
            coefficient = Method('CoefficientType< i > &', 'coefficient', targs=['std::size_t i'], code=return_('std::get< i >( coefficients_ )'))
        else:
            coefficient = Method('CoefficientType< i > &', 'coefficient', targs=['std::size_t i', 'Side side = Side::in'], code=return_('std::get< i >( coefficients_[ side ] )'))
        result += [coefficient.variant('const CoefficientType< i > &coefficient', const=True), coefficient]

        result.append(AccessModifier('private'))

        initCoefficients = Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['const EntityType &entity', 'std::tuple< Coefficients... > &coefficients', 'std::index_sequence< i... >'], static=True)
        initCoefficients.append('std::ignore = std::make_tuple( (std::get< i >( coefficients ).init( entity ), i)... );')
        result.append(initCoefficients)

        constructConstants = Method('void', 'constructConstants', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        constructConstants.append('std::ignore = std::make_tuple( (std::get< i >( constants_ ) = std::make_shared< ConstantsType< i > >(), i)... );')
        result.append(constructConstants)

        if self.skeleton is None:
            result.append(Declaration(Variable('EntityType', 'entity_')))
            result.append(Declaration(Variable('std::tuple< Coefficients... >', 'coefficients_')))
        else:
            result.append(Declaration(Variable('std::array< EntityType, 2 >', 'entity_')))
            result.append(Declaration(Variable('std::array< std::tuple< Coefficients... >, 2 >', 'coefficients_')))
        result.append(Declaration(Variable('ConstantsTupleType', 'constants_')))
        if self.vars is not None:
            result += self.vars

        return result

    def write(self, sourceWriter, name='Integrands', targs=[]):
        code = Struct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']))
        code.append(self.pre())
        code.append(self.main())
        code.append(self.post())
        sourceWriter.emit(code)


def generateCode(predefined, testFunctions, tensorMap, coefficients, tempVars=True):
    # build list of all expressions to compile
    expressions = []
    for phi in testFunctions:
        if (phi,) not in tensorMap:
            continue
        tensor = tensorMap[(phi,)]
        keys = tensor.keys()
        expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars)

    # extract generated code for expressions and build values
    values = []
    for phi in testFunctions:
        value = tensors.fill(phi.ufl_shape, '0')
        if (phi,) in tensorMap:
            tensor = tensorMap[(phi,)]
            keys = tensor.keys()
            for i, r in zip(keys, results[:len(keys)]):
                value = tensors.setItem(value, i, r)
            results = results[len(keys):]
        values += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values


def generateLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, coefficients, tempVars=True):
    # build list of all expressions to compile
    expressions = []
    for phi in testFunctions:
        for psi in trialFunctions:
            if (phi, psi) not in tensorMap:
                continue
            tensor = tensorMap[(phi, psi)]
            keys = tensor.keys()
            expressions += [tensor[i] for i in keys]

    # compile all expressions at once
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars)

    # extract generated code for expressions and build values
    values = []
    for phi in testFunctions:
        value = tensors.fill(phi.ufl_shape, None)
        for idx in range(len(trialFunctions)):
            psi = trialFunctions[idx]
            varpsi = Variable('std::tuple< RangeType, JacobianRangeType >', 'phi')[idx]
            if (phi, psi) in tensorMap:
                tensor = tensorMap[(phi, psi)]
                keys = tensor.keys()
                for ij, r in zip(keys, results[:len(keys)]):
                    if isinstance(tensor[ij], Zero):
                        continue
                    i = ij[:len(phi.ufl_shape)]
                    j = ij[len(phi.ufl_shape):]
                    if isinstance(tensor[ij], IntValue) and int(tensor[ij]) == 1:
                        r = varpsi[j]
                    else:
                        r = r * varpsi[j]
                    s = tensors.getItem(value, i)
                    s = r if s is None else s + r
                    value = tensors.setItem(value, i, s)
                results = results[len(keys):]
        value = tensors.apply(lambda v : makeExpression(0) if v is None else v, phi.ufl_shape, value)
        values += [tensors.reformat(lambda row: InitializerList(*row), phi.ufl_shape, value)]

    return preamble, values


def generateUnaryCode(predefined, testFunctions, tensorMap, coefficients, tempVars=True):
    preamble, values = generateCode(predefined, testFunctions, tensorMap, coefficients, tempVars)
    return preamble + [return_(construct('ValueType', *values))]


def generateUnaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, coefficients, tempVars=True):
    if tensorMap is None:
        return [return_(lambda_(args=['const ValueType &phi'], code=return_(construct('ValueType', 0, 0))))]

    preamble, values = generateLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, coefficients, tempVars)
    capture = extractVariables(values) - {Variable('std::tuple< RangeType, JacobianRangeType >', 'phi')}
    return preamble + [return_(lambda_(capture=capture, args=['const ValueType &phi'], code=return_(construct('ValueType', *values))))]


def generateBinaryCode(predefined, testFunctions, tensorMap, coefficients, tempVars=True):
    restrictedTestFunctions = [phi('-') for phi in testFunctions] + [phi('+') for phi in testFunctions]
    preamble, values = generateCode(predefined, restrictedTestFunctions, tensorMap, coefficients, tempVars=True)
    return preamble + [return_(make_pair(construct('ValueType', *values[:len(testFunctions)]), construct('ValueType', *values[len(testFunctions):])))]


def generateBinaryLinearizedCode(predefined, testFunctions, trialFunctions, tensorMap, coefficients, tempVars=True):
    restrictedTestFunctions = [phi('-') for phi in testFunctions] + [phi('+') for phi in testFunctions]
    restrictedTrialFunctions = [psi('-') for psi in trialFunctions] + [psi('+') for psi in trialFunctions]
    if tensorMap is None:
        return [return_(lambda_(args=['const ValueType &phi'], code=return_(construct('ValueType', 0, 0))))]

    preamble, values = generateLinearizedCode(predefined, restrictedTestFunctions, restrictedTrialFunctions, tensorMap, coefficients, tempVars)
    tensorIn = lambda_(args=['const ValueType &phiIn'], code=return_(construct('ValueType', *values[:len(testFunctions)])))
    tensorOut = lambda_(args=['const ValueType &phiOut'], code=return_(construct('ValueType', *values[len(testFunctions):])))
    return preamble + [return_(make_pair(tensorIn, tensorOut))]

def compileUFL(equation, tempVars=True):
    form = equation.lhs - equation.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    phi = form.arguments()[0]
    dphi = Grad(phi)
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = Grad(u)
    d2u = Grad(du)
    ubar = Coefficient(u.ufl_function_space())
    dubar = Grad(ubar)
    d2ubar = Grad(dubar)

    integrands = Integrands(dimRange, form.signature())
    try:
        integrands.field = u.ufl_function_space().ufl_element().field()
    except AttributeError:
        pass

    idxConst = 0
    idxCoeff = 0
    for coefficient in set(form.coefficients()):
        try:
            name = coefficient.str()
        except:
            name = str(coefficient)
        if coefficient.is_cellwise_constant():
            dimRange = 1 if coefficient.ufl_shape==() else coefficient.ufl_shape[0]
            integrands.coefficients.append({'name': name, 'number' : idxConst, 'counter' : coefficient.count(), 'dimRange' : dimRange, 'constant': True, 'field': None})
            idxConst += 1
        else:
            field = coefficient.ufl_function_space().ufl_element().field()
            integrands.coefficients.append({'name': name, 'number' : idxCoeff, 'counter' : coefficient.count(), 'dimRange' : coefficient.ufl_shape[0], 'constant': False, 'field': field})
            idxCoeff += 1

    integrals = splitForm(form, [phi])

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))
    linearizedIntegrals = splitForm(dform, [phi, u])

    if not set(integrals.keys()) <= {'cell', 'exterior_facet', 'interior_facet'}:
        raise Exception('unknown integral encountered in ' + str(set(integrals.keys())) + '.')

    if 'cell' in integrals.keys():
        arg = Variable('std::tuple< RangeType, JacobianRangeType >', 'u')
        predefined = {u: arg[0], du: arg[1]}
        integrands.interior = generateUnaryCode(predefined, [phi, dphi], integrals['cell'], integrands.coefficients, tempVars)
        predefined = {ubar: arg[0], dubar: arg[1]}
        integrands.linearizedInterior = generateUnaryLinearizedCode(predefined, [phi, dphi], [u, du], linearizedIntegrals.get('cell'), integrands.coefficients, tempVars)

    if 'exterior_facet' in integrals.keys():
        predefined = {u: arg[0], du: arg[1]}
        integrands.boundary = generateUnaryCode(predefined, [phi, dphi], integrals['exterior_facet'], integrands.coefficients, tempVars);
        predefined = {ubar: arg[0], dubar: arg[1]}
        integrands.linearizedBoundary = generateUnaryLinearizedCode(predefined, [phi, dphi], [u, du], linearizedIntegrals.get('exterior_facet'), integrands.coefficients, tempVars)

    if 'interior_facet' in integrals.keys():
        argIn = Variable('std::tuple< RangeType, JacobianRangeType >', 'uIn')
        argOut = Variable('std::tuple< RangeType, JacobianRangeType >', 'uOut')
        predefined = {u('-'): argIn[0], du('-'): argIn[1], u('+'): argOut[0], du('+'): argOut[1]}
        integrands.skeleton = generateBinaryCode(predefined, [phi, dphi], integrals['interior_facet'], integrands.coefficients, tempVars)
        predefined = {ubar('-'): argIn[0], dubar('-'): argIn[1], ubar('+'): argOut[0], dubar('+'): argOut[1]}
        integrands.linearizedSkeleton = generateBinaryLinearizedCode(predefined, [phi, dphi], [u, du], linearizedIntegrals.get('interior_facet'), integrands.coefficients, tempVars)

    return integrands
