from __future__ import division, print_function

from ufl import Coefficient, Form
from ufl import action, derivative
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.differentiation import Grad

from dune.source.cplusplus import AccessModifier, Constructor, EnumClass, Method, Struct, TypeAlias, Variable
from dune.source.cplusplus import SourceWriter

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

        result.append(Variable("const int dimRange", str(self.dimRange), static=True))
        result.append(Variable("const int dimDomain", "GridPartType::dimensionworld", static=True))
        result.append(Variable("const int dimLocal", "GridPartType::dimension", static=True))

        result.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        result.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))
        result.append(TypeAlias("FunctionSpaceType", "Dune::Fem::FunctionSpace< double, RangeFieldType, dimDomain, dimRange >"))

        for s in ["DomainType", "RangeType", "JacobianRangeType", "HessianRangeType"]:
            result.append(TypeAlias(s, "typename FunctionSpaceType::" + s))

        if self.coefficients:
            constants = [("std::shared_ptr< Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " > >") for c in self.coefficients if c['constant']]
            result.append(TypeAlias("ConstantsTupleType", "std::tuple< " + ", ".join(constants) + " >"))
            coefficients = [('Dune::Fem::FunctionSpace< double,' + SourceWriter.cpp_fields(c['field']) + ', dimDomain, ' + str(c['dimRange']) + ' >') for c in self.coefficients if not c['constant']]
            result.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " +", ".join(coefficients) + " >"))

            result.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantsTupleType >::element_type", targs=["std::size_t i"]))
            result.append(Variable("const std::size_t numCoefficients", "std::tuple_size< CoefficientFunctionSpaceTupleType >::value", static=True))
            result.append(TypeAlias("CoefficientFunctionSpaceType", "typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type", targs=["std::size_t i"]))
            for s in ["RangeType", "JacobianRangeType", "HessianRangeType"]:
                result.append(TypeAlias("Coefficient" + s, "typename CoefficientFunctionSpaceType< i >::" + s, targs=["std::size_t i"]))
        else:
            result.append(Variable("const std::size_t numCoefficients", "0u", static=True))
            result.append(TypeAlias("ConstantsTupleType", "std::tuple<>"))

        result.append(TypeAlias('CoefficientType', 'typename std::tuple_element< i, std::tuple< Coefficients... > >::type', targs=['std::size_t i']))
        result.append(TypeAlias('ConstantsType', 'typename std::tuple_element< i, ConstantsTupleType >::type::element_type', targs=['std::size_t i']))

        if self.skeleton is not None:
            result.append(EnumClass('Side', ['in = 0u', 'out = 1u'], 'std::size_t'))
            inside = '[ Side::in ]'
        else:
            inside = ''

        result.append(Constructor(code="constructConstants( std::make_index_sequence< std::tuple_size< ConstantsTupleType >::value >() );"))

        initEntity = Method('bool init', args=['const EntityType &entity'], const=True)
        initEntity.append('entity_' + inside + ' = entity;')
        initEntity.append('initCoefficients( entity_' + inside + ', coefficients_' + inside + ', std::make_index_sequence< numCoefficients >() );')
        initEntity.append(self.init)
        initEntity.append('return true;')
        result.append(initEntity)

        initIntersection = Method('bool init', args=['const IntersectionType &intersection'], const=True)
        if self.skeleton is None:
            initIntersection.append('return (intersection.boundary() && init( intersection.inside() ));')
        else:
            initIntersection.append('entity_[ Side::in ] = intersection.inside();')
            initIntersection.append('initCoefficients( entity_[ Side::in ], coefficients_[ Side::in ], std::make_index_sequence< numCoefficients >() );')
            initIntersection.append('if( intersection.neighbor() )')
            initIntersection.append('{')
            initIntersection.append('  entity_[ Side::out ] = intersection.outside();')
            initIntersection.append('  initCoefficients( entity_[ Side::out ], coefficients_[ Side::out ], std::make_index_sequence< numCoefficients >() );')
            initIntersection.append('}')
            initIntersection.append('return true;')
        result.append(initIntersection)

        return result

    def main(self):
        result = []

        if self.interior is not None:
            result.append(Method('ValueType interior', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.interior, const=True))
            #result.append(Method('auto linearizedInterior', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.linearizedInterior, const=True))

        if self.boundary is not None:
            result.append(Method('ValueType boundary', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.boundary, const=True))
            #result.append(Method('auto linearizedBoundary', targs=['class Point'], args=['const Point &x', 'const ValueType &u'], code=self.linearizedBoundary, const=True))

        if self.skeleton is not None:
            result.append(Method('std::pair< ValueType, ValueType > skeleton', targs=['class Point'], args=['const Point &xIn', 'const ValueType &uIn', 'const Point &xOut', 'const ValueType &uOut'], code=self.skeleton, const=True))
            #result.append(Method('auto linearizedSkeleton', targs['class Point'], args=['const Point &xIn', 'const ValueType &uIn', 'const Point &xOut', 'const ValueType &uOut'], code=self.linearizedSkeleton, const=True))

        return result

    def post(self):
        result = []

        constant = Method('ConstantsType< i > &constant', targs=['std::size_t i'], code='return *std::get< i >( constants_ );')
        result += [constant.variant('const ConstantsType< i > &constant', const=True), constant]

        if self.skeleton is None:
            coefficient = Method('CoefficientType< i > &coefficient', targs=['std::size_t i'], code='return std::get< i >( coefficients_ );')
        else:
            coefficient = Method('CoefficientType< i > &coefficient', targs=['std::size_t i', 'Side side = Side::in'], code='return std::get< i >( coefficients_[ side ] );')
        result += [coefficient.variant('const CoefficientType< i > &coefficient', const=True), coefficient]

        result.append(AccessModifier('private'))

        initCoefficients = Method('void initCoefficients', targs=['std::size_t... i'], args=['const EntityType &entity', 'std::tuple< Coefficients... > &coefficients', 'std::index_sequence< i... >'], static=True)
        initCoefficients.append('std::ignore = std::make_tuple( (std::get< i >( coefficients ).init( entity ), i)... );')
        result.append(initCoefficients)

        constructConstants = Method('void constructConstants', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        constructConstants.append('std::ignore = std::make_tuple( (std::get< i >( constants_ ) = std::make_shared< ConstantsType< i > >(), i)... );')
        result.append(constructConstants)

        if self.skeleton is None:
            result.append(Variable('EntityType entity_'))
            result.append(Variable('std::tuple< Coefficients... > coefficients_'))
        else:
            result.append(Variable('std::array< EntityType, 2 > entity_'))
            result.append(Variable('std::array< std::tuple< Coefficients... >, 2 > coefficients_'))
        result.append(Variable('ConstantsTupleType constants_'))
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
        values += [tensors.reformat(lambda row: '{ ' + ', '.join(row) + ' }', phi.ufl_shape, value)]

    return preamble, values


def generateUnaryCode(predefined, testFunctions, tensorMap, coefficients, tempVars=True):
    preamble, values = generateCode(predefined, testFunctions, tensorMap, coefficients, tempVars)
    return preamble + ['return ValueType( ' + ', '.join(values) + ' );']


def generateBinaryCode(predefined, testFunctions, tensorMap, coefficients, tempVars=True):
    restrictedTestFunctions = [phi('-') for phi in testFunctions] + [phi('+') for phi in testFunctions]
    preamble, values = generateCode(predefined, restrictedTestFunctions, tensorMap, coefficients, tempVars=True)
    return preamble + ['return std::make_pair( ValueType( ' + ', '.join(values) + ' ), ValueType( ' + ', '.join(values) + ' ) );']


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

    #if integrals.keys() != linearizedIntegrals.keys():
    #    raise Exception('linearization contains different integrals.')
    if not set(integrals.keys()) <= {'cell', 'exterior_facet', 'interior_facet'}:
        raise Exception('unknown integral encountered in ' + str(set(integrals.keys())) + '.')

    if 'cell' in integrals.keys():
        predefined = {u: 'std::get< 0 >( u )', du: 'std::get< 1 >( u )'}
        integrands.interior = generateUnaryCode(predefined, [phi, dphi], integrals['cell'], integrands.coefficients, tempVars);

    if 'exterior_facet' in integrals.keys():
        predefined = {u: 'std::get< 0 >( u )', du: 'std::get< 1 >( u )'}
        integrands.boundary = generateUnaryCode(predefined, [phi, dphi], integrals['exterior_facet'], integrands.coefficients, tempVars);

    if 'interior_facet' in integrals.keys():
        for key, value in integrals['interior_facet'].items():
            print('----- ' + ' | '.join([str(k) for k in key]) + ' -----')
            for k in value.keys():
                print(', '.join([str(i) for i in k]) + ': ' + str(value[k]))
        predefined = {u('-'): 'std::get< 0 >( uIn )', du('-'): 'std::get< 1 >( uIn )', u('+'): 'std::get< 0 >( uOut )', du('+'): 'std::get< 0 >( uOut )'}
        integrands.skeleton = generateBinaryCode(predefined, [phi, dphi], integrals['interior_facet'], integrands.coefficients, tempVars)

    return integrands
