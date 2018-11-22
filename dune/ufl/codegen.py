from __future__ import division, print_function

import re

from ufl.algorithms import expand_indices
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.corealg.map_dag import map_expr_dags
from ufl.corealg.multifunction import MultiFunction
from ufl.argument import Argument
from ufl.coefficient import Coefficient
from ufl.differentiation import Grad
from ufl.core.multiindex import FixedIndex, MultiIndex

from dune.source.builtin import get, hybridForEach, make_pair, make_index_sequence, make_shared
import dune.source.cplusplus as cplusplus
from dune.source.cplusplus import ConditionalExpression, Declaration, Using, Variable

from dune.source.cplusplus import AccessModifier, Declaration, Constructor, EnumClass, Include, InitializerList, Method, Struct, TypeAlias, UnformattedExpression, Variable
from dune.source.cplusplus import assign, construct, coordinate, dereference, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_

from .applyrestrictions import applyRestrictions

def translateIndex(index):
    if isinstance(index, (tuple, MultiIndex)):
        return ''.join([translateIndex(i) for i in index])
    elif isinstance(index, (int, FixedIndex)):
        return '[ ' + str(index) + ' ]'
    else:
        raise Exception('Index type not supported: ' + repr(index))


class CodeGenerator(MultiFunction):
    def __init__(self, predefined, coefficients, tempVars):
        MultiFunction.__init__(self)
        self.using = set()
        self.predefined = predefined
        self.coefficients = [] if coefficients is None else coefficients
        self.code = []
        self.tempVars = tempVars

    def _require_predefined(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            raise Exception('%s not available for this expression.' % expr._ufl_class_.__name__)

    def abs(self, expr, x):
        self.using.add(Using(cplusplus.abs_))
        return self._makeTmp(cplusplus.abs_(x))

    def AndCondition(self, expr, left, right):
        self.using.add(Using(cplusplus.and_))
        return self._makeTmp(cplusplus.and_(left, right))

    def argument(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            raise Exception('Unknown argument: ' + str(expr.number()))

    def atan(self, expr, x):
        self.using.add(Using(cplusplus.atan))
        return self._makeTmp(cplusplus.atan(x))

    def atan_2(self, expr, x, y):
        self.using.add(Using(cplusplus.atan2))
        return self._makeTmp(cplusplus.atan2(x, y))

    cell_volume = _require_predefined

    def coefficient(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        # print('Warning: ' + ('Constant ' if expr.is_cellwise_constant() else 'Coefficient ') + str(expr) + ' not predefined.')
        idx = str(self._getNumber(expr))
        if expr.is_cellwise_constant():
            var = Variable('const ConstantsRangeType< ' + idx + ' >', 'cc' + idx)
            self.code.append(Declaration(var, 'constant< ' + idx + ' >()'))
        else:
            var = Variable('CoefficientRangeType< ' + idx + ' >', 'c' + idx)
            self.code.append(Declaration(var))
            self.code.append('coefficient< ' + idx + ' >().evaluate( x, c' + idx + ' );')
        return var

    def conditional(self, expr, cond, true, false):
        return ConditionalExpression('auto', cond, true, false)

    def cos(self, expr, x):
        self.using.add(Using(cplusplus.cos))
        return self._makeTmp(cplusplus.cos(x))

    def cosh(self, expr, x):
        self.using.add(Using(cplusplus.cosh))
        return self._makeTmp(cplusplus.cosh(x))

    def division(self, expr, x, y):
        return self._makeTmp(x / y)

    def eq(self, expr, left, right):
        return left == right

    def exp(self, expr, x):
        self.using.add(Using(cplusplus.exp))
        return self._makeTmp(cplusplus.exp(x))

    facet_area = _require_predefined
    facet_normal = _require_predefined

    def float_value(self, expr):
        return cplusplus.makeExpression(expr.value())

    def ge(self, expr, left, right):
        return left >= right

    def gt(self, expr, left, right):
        return left > right

    def grad(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        operand = expr.ufl_operands[0]
        if isinstance(operand, Coefficient):
            idx = str(self._getNumber(operand))
            var = Variable('CoefficientJacobianRangeType< ' + idx + ' >', 'dc' + idx)
            self.code.append(Declaration(var))
            self.code.append('coefficient< ' + idx + ' >().jacobian( x, ' + var.name + ' );')
            return var
        elif isinstance(operand, Grad):
            operand = operand.ufl_operands[0]
            if isinstance(operand, Coefficient):
                idx = str(self._getNumber(operand))
                var = Variable('CoefficientHessianRangeType< ' + idx + ' >', 'd2c' + idx)
                self.code.append(Declaration(var))
                self.code.append('coefficient< ' + idx + ' >().hessian( x, ' + var.name + ' );')
                return var
            elif isinstance(operand, Argument):
                raise Exception('Unknown argument: ' + str(operand))
            else:
                raise Exception('CodeGenerator does not allow for third order derivatives, yet.')
        elif isinstance(operand, Argument):
            raise Exception('Unknown argument: ' + str(operand))
        else:
            raise Exception('Cannot compute gradient of ' + repr(expr))

    def indexed(self, expr, operand, index):
        for i in index:
            operand = operand[int(i)]
        return operand

    def le(self, expr, left, right):
        return left <= right

    def ln(self, expr, x):
        self.using.add(Using(cplusplus.log))
        return self._makeTmp(cplusplus.log(x))

    def lt(self, expr, left, right):
        return left < right

    max_cell_edge_length = _require_predefined
    max_facet_edge_length = _require_predefined

    def MaxValue(self, expr, left, right):
        self.using.add(Using(cplusplus.max_))
        return self._makeTmp(cplusplus.max_(left, right))
    def max_value(self, expr, left, right):
        self.using.add(Using(cplusplus.max_))
        return self._makeTmp(cplusplus.max_(left, right))

    min_cell_edge_length = _require_predefined
    min_facet_edge_length = _require_predefined

    def MinValue(self, expr, left, right):
        self.using.add(Using(cplusplus.min_))
        return self._makeTmp(cplusplus.min_(left, right))
    def min_value(self, expr, left, right):
        self.using.add(Using(cplusplus.min_))
        return self._makeTmp(cplusplus.min_(left, right))

    def multi_index(self, expr):
        return expr

    int_value = float_value

    def restricted(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        operand = expr.ufl_operands[0]
        if isinstance(operand, Coefficient) and operand.ufl_element().family() == "Real":
            return self.coefficient(operand)
        raise Exception('Cannot compute restriction of ' + str(operand))

    def product(self, expr, x, y):
        return self._makeTmp(x * y)

    def power(self, expr, x, y):
        self.using.add(Using(cplusplus.pow_))
        return self._makeTmp(cplusplus.pow_(x, y))

    def sin(self, expr, x):
        self.using.add(Using(cplusplus.sin))
        return self._makeTmp(cplusplus.sin(x))

    def sinh(self, expr, x):
        self.using.add(Using(cplusplus.sinh))
        return self._makeTmp(cplusplus.sinh(x))

    def sqrt(self, expr, x):
        self.using.add(Using(cplusplus.sqrt))
        return self._makeTmp(cplusplus.sqrt(x))

    def spatial_coordinate(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            self.using.add(Using(cplusplus.coordinate))
            var = Variable('const auto', 'y')
            self.code.append(Declaration(var, 'entity().geometry().global( coordinate( x ) )'))
        return var

    def sum(self, expr, x, y):
        return self._makeTmp(x + y)

    def tan(self, expr, x):
        self.using.add(Using(cplusplus.tan))
        return self._makeTmp(cplusplus.tan(x))

    def tanh(self, expr, x):
        self.using.add(Using(cplusplus.tanh))
        return self._makeTmp(cplusplus.tanh(x))

    def zero(self, expr):
        return cplusplus.makeExpression(0)

    def _getNumber(self, expr):
        try:
            name = getattr(expr, "name")
        except AttributeError:
            name = str(expr)
        e = [ee for ee in self.coefficients if ee["name"] == name]
        if len(e) > 1:
            raise KeyError('two coefficients provided with same name')
        if len(e) == 0:
            raise KeyError('coefficient provided with no name')
        return e[0]["number"]

    def _makeTmp(self, cexpr, tempVars=None):
        if isinstance(cexpr, Variable):
            return cexpr
        if tempVars is None:
            tempVars = self.tempVars
        if tempVars:
            cppType = None
            if isinstance(cexpr, cplusplus.Expression):
                cppType = cexpr.cppType
            if cppType is None:
                cppType = 'const auto'
            var = Variable(cppType, 'tmp' + str(len(self.code)))
            self.code.append(Declaration(var, cexpr))
            return var
        else:
            return cexpr


def generateCode(predefined, expressions, coefficients=None, tempVars=True):
    expressions = [applyRestrictions(expand_indices(apply_derivatives(apply_algebra_lowering(expr)))) for expr in expressions]
    generator = CodeGenerator(predefined, coefficients, tempVars)
    results = map_expr_dags(generator, expressions)
    l = list(generator.using)
    l.sort()
    return l + generator.code, results


######################################################################

def fieldVectorType(shape, field = None, useScalar = False):
    if isinstance(shape, Coefficient):
        if field is not None:
            raise ValueError("Cannot specify field type for coefficients")
        try:
            field = shape.ufl_function_space().field()
        except AttributeError:
            field = 'double'
        shape = shape.ufl_shape
    else:
        field = 'double' if field is None else field

    field = 'std::complex< double >' if field == 'complex' else field

    if not isinstance(shape, tuple):
        raise ValueError("Shape must be a tuple.")
    dimRange = (1 if len(shape) == 0 else shape[0])

    if dimRange == 1 and useScalar:
        return field
    else:
        return 'Dune::FieldVector< ' + field + ', ' + str(dimRange) + ' >'

class ModelClass():
    def __init__(self, uflExpr):
        coefficients = set()
        for expr in uflExpr:
            try:
                coefficients |= set(expr.coefficients())
            except:
                _, cc = extract_arguments_and_coefficients(expr)
                coefficients |= set(cc)

        self.constantList = [c for c in coefficients if c.is_cellwise_constant()]
        self.coefficientList = sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())

        constants=(fieldVectorType(c,useScalar=True) for c in self.constantList)
        coefficients=(fieldVectorType(c) for c in self.coefficientList)
        constantNames=(getattr(c, 'name', None) for c in self.constantList)
        coefficientNames=(getattr(c, 'name', None) for c in self.coefficientList)
        parameterNames=(getattr(c, 'parameter', None) for c in self.constantList)

        self._constants = [] if constants is None else list(constants)
        self._coefficients = [] if coefficients is None else list(coefficients)

        self._constantNames = [None,] * len(self._constants) if constantNames is None else list(constantNames)
        self._constantNames = ['constant' + str(i) if n is None else n for i, n in enumerate(self._constantNames)]
        if len(self._constantNames) != len(self._constants):
            raise ValueError("Length of constantNames must match length of constants")
        invalidConstants = [n for n in self._constantNames if n is not None and re.match('^[a-zA-Z_][a-zA-Z0-9_]*$', n) is None]
        if invalidConstants:
            raise ValueError('Constant names are not valid C++ identifiers:' + ', '.join(invalidCoefficients) + '.')

        self._coefficientNames = [None,] * len(self._coefficients) if coefficientNames is None else list(coefficientNames)
        self._coefficientNames = ['coefficient' + str(i) if n is None else n for i, n in enumerate(self._coefficientNames)]
        if len(self._coefficientNames) != len(self._coefficients):
            raise ValueError("Length of coefficientNames must match length of coefficients")
        invalidCoefficients = [n for n in self._coefficientNames if n is not None and re.match('^[a-zA-Z_][a-zA-Z0-9_]*$', n) is None]
        if invalidCoefficients:
            raise ValueError('Coefficient names are not valid C++ identifiers:' + ', '.join(invalidCoefficients) + '.')

        self._parameterNames = [None,] * len(self._constants) if parameterNames is None else list(parameterNames)
        if len(self._parameterNames) != len(self._constants):
            raise ValueError("Length of parameterNames must match length of constants")
        self._derivatives = [('RangeType', 'evaluate'), ('JacobianRangeType', 'jacobian'), ('HessianRangeType', 'hessian')]

    @property
    def constantTypes(self):
        return [n[0].upper() + n[1:] for n in self._constantNames]

    @property
    def constantNames(self):
        return [n[0].lower() + n[1:] for n in self._constantNames]

    @property
    def coefficientTypes(self):
        return [n[0].upper() + n[1:] for n in self._coefficientNames]

    @property
    def coefficientNames(self):
        return [n[0].lower() + n[1:] for n in self._coefficientNames]

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def coefficient(self, idx, x, side=None):
        targs = [str(idx)]
        if side is not None:
            targs.append(side)
        return (UnformattedExpression('typename CoefficientFunctionSpaceType< ' + str(idx) + ' >::' + t, n + 'Coefficient< ' + ', '.join(targs) + ' >( ' + x + ' )') for t, n in self._derivatives)

    def _predefineCoefficients(self, predefined, x, side=None):
        for idx, coefficient in enumerate(self.coefficientList):
            for derivative in self.coefficient(idx, x, side=side):
                if side is None:
                    predefined[coefficient] = derivative
                elif side == 'Side::in':
                    predefined[coefficient('+')] = derivative
                elif side == 'Side::out':
                    predefined[coefficient('-')] = derivative
                coefficient = Grad(coefficient)
    def predefineCoefficients(self, predefined, skeleton=False):
        predefined.update({c: self.constant(i) for i, c in enumerate(self.constantList)})
        if skeleton is False:
            self._predefineCoefficients(predefined, 'x')
        else:
            self._predefineCoefficients(predefined, 'xIn', 'Side::in')
            self._predefineCoefficients(predefined, 'xOut', 'Side::out')



# model.name = "Integrands"
# model.targs = ['class GridPart']
# model.gridPartType = TypeAlias("GridPartType", "GridPart")
# model.ctor_args = []
# model.ctor_init = []
# model.skeleton (None is no intersecton terms are needed)
# model.init (code to be added to init(Entity) method
# model.vars (further class variables)

####

# model._coffieicnts      # c++ type of coefficient range
# model.coefficientTypes  # type of coefficient (i.e. template argument)
# model.coefficientNames  # name of variable passed to ctor

# model._constants      # C++ types
# model.constantTypes   # type alias
# model.constantNames   # names

# model.parameterNames

# model._derivatives

def generateModelClass(model):
    code = Struct(model.name, targs=(model.targs + ['class ' + n for n in model.coefficientTypes]))

    code.append(model.gridPartType)

    code.append(TypeAlias("EntityType", "typename GridPartType::template Codim< 0 >::EntityType"))
    code.append(TypeAlias("IntersectionType", "typename GridPartType::IntersectionType"))

    code.append(TypeAlias("GlobalCoordinateType", "typename EntityType::Geometry::GlobalCoordinate"))

    for type, alias in zip(model._constants, model.constantTypes):
        code.append(TypeAlias(alias, type))
    constants = ["std::shared_ptr< " + c + " >" for c in model.constantTypes]
    if constants:
        code.append(TypeAlias("ConstantTupleType", "std::tuple< " + ", ".join(constants) + " >"))
        code.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantTupleType >::element_type", targs=["std::size_t i"]))
    else:
        code.append(TypeAlias("ConstantTupleType", "std::tuple<>"))

    if model._coefficients:
        coefficientSpaces = [('Dune::Fem::GridFunctionSpace< GridPartType, ' + c + ' >') for c in model._coefficients]
        code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " +", ".join(coefficientSpaces) + " >"))
        code.append(TypeAlias('CoefficientTupleType', 'std::tuple< ' + ', '.join(model.coefficientTypes) + ' >'))

        code.append(TypeAlias("CoefficientFunctionSpaceType", "std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >", targs=["std::size_t i"]))
        for s in ["RangeType", "JacobianRangeType"]:
            code.append(TypeAlias("Coefficient" + s, "typename CoefficientFunctionSpaceType< i >::" + s, targs=["std::size_t i"]))
    else:
        code.append(TypeAlias("CoefficientTupleType", "std::tuple<>"))

    code.append(TypeAlias('CoefficientType', 'std::tuple_element_t< i, CoefficientTupleType >', targs=['std::size_t i']))
    code.append(TypeAlias('ConstantType', 'typename std::tuple_element_t< i, ConstantTupleType >::element_type', targs=['std::size_t i']))

    if model.skeleton is not None:
        code.append(EnumClass('Side', ['in = 0u', 'out = 1u'], 'std::size_t'))
        inside = '[ static_cast< std::size_t >( Side::in ) ]'
    else:
        inside = ''

    if model.skeleton is None:
        entity_ = Variable('EntityType', 'entity_')
        insideEntity = entity_
    else:
        entity_ = Variable('std::array< EntityType, 2 >', 'entity_')
        insideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::in )')]
        outsideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::out )')]
    intersection_ = Variable('IntersectionType', 'intersection_')

    constants_ = Variable("ConstantTupleType", "constants_")
    coefficientsTupleType = 'std::tuple< ' + ', '.join('Dune::Fem::ConstLocalFunction< ' + n + ' >' for n in model.coefficientTypes) + ' >'
    if model.skeleton is None:
        coefficients_ = Variable(coefficientsTupleType, 'coefficients_')
    else:
        coefficients_ = Variable('std::array< ' + coefficientsTupleType + ', 2 >', 'coefficients_')

    # generate code for constructor
    arg_param = Variable('const Dune::Fem::ParameterReader &', 'parameter')
    args = model.ctor_args + [Variable('const ' + t + ' &', n) for t, n in zip(model.coefficientTypes, model.coefficientNames)]
    if model._coefficients:
        init = ['Dune::Fem::ConstLocalFunction< ' + n + ' >( ' + p + ' )' for n, p in zip(model.coefficientTypes, model.coefficientNames)]
        if model.skeleton is None:
            init = ["coefficients_( " + ", ".join(init) + " )"]
        else:
            init = ['coefficients_{{ ' + coefficientsTupleType + '( ' + ', '.join(init) + ' ), ' + coefficientsTupleType + '( ' + ', '.join(init) + ' ) }}']
    else:
        init = []
    init = model.ctor_init + init
    args.append(Declaration(arg_param, initializer=UnformattedExpression('const ParameterReader &', 'Dune::Fem::Parameter::container()')))
    constructor = Constructor(args=args, init=init)
    for idx, cppType in enumerate(model.constantTypes):
        constructor.append(assign(get(idx)(constants_), make_shared(cppType)()))
    for idx, (name, cppType) in enumerate(zip(model._parameterNames, model.constantTypes)):
        if name is not None:
            constructor.append(assign(dereference(get(idx)(constants_)), UnformattedExpression('auto', arg_param.name + '.getValue< ' + cppType + ' >( "' + name + '" )', uses=[arg_param])))
    code.append(constructor)

    entity = Variable('const EntityType &', 'entity')
    initEntity = Method('bool', 'init', args=[entity])
    initEntity.append(assign(insideEntity, entity))
    if model.skeleton is None:
        for i, c in enumerate(model._coefficients):
            initEntity.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + ' ).init( entity )', uses=[entity, coefficients_]))
    else:
        for i, c in enumerate(model._coefficients):
            initEntity.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::in ) ] ).init( entity )', uses=[entity, coefficients_]))
    initEntity.append(model.init)
    initEntity.append(return_(True))
    code.append(initEntity)

    intersection = Variable('const IntersectionType &', 'intersection')
    initIntersection = Method('bool', 'init', args=[intersection])
    initIntersection.append(assign(intersection_, intersection))
    if model.skeleton is None:
        initIntersection.append(return_('(intersection.boundary() && init( intersection.inside() ))'))
    else:
        initIntersection.append(assign(insideEntity, UnformattedExpression('EntityType', 'intersection.inside()')))
        for i, c in enumerate(model._coefficients):
            initIntersection.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::in ) ] ).init( entity_[ static_cast< std::size_t >( Side::in ) ] )', uses=[coefficients_]))
        initIntersection.append('if( intersection.neighbor() )')
        initIntersection.append('{')
        initIntersection.append('  entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();')
        for i, c in enumerate(model._coefficients):
            initIntersection.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::out ) ] ).init( entity_[ static_cast< std::size_t >( Side::out ) ] )', uses=[coefficients_]))
        initIntersection.append('}')
        initIntersection.append(return_(True))
    code.append(initIntersection)

################################
    model.methods(code)
################################

    code.append(Method('const ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_))), const=True))
    code.append(Method('ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_)))))

    for i, (t, n) in enumerate(zip(model.constantTypes, model.constantNames)):
        code.append(Method('const ' + t + ' &', n, code=return_(dereference(get(i)(constants_))), const=True))
        code.append(Method(t + ' &', n, code=return_(dereference(get(i)(constants_)))))

    code.append(Method('const EntityType &', 'entity', const=True, code=return_(insideEntity)))

    code.append(AccessModifier('private'))

    if model._coefficients:
        for cppType, name in model._derivatives:
            var = Variable('typename CoefficientFunctionSpaceType< i >::' + cppType, 'result')
            if model.skeleton is None:
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

    if model.vars is not None:
        code += model.vars

    return code
