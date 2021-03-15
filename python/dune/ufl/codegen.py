from __future__ import division, print_function

import re

from ufl import replace
from ufl.algorithms import expand_indices
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.algorithms.apply_algebra_lowering import apply_algebra_lowering
from ufl.corealg.map_dag import map_expr_dags
from ufl.corealg.multifunction import MultiFunction
from ufl.argument import Argument
from ufl.coefficient import Coefficient
from ufl.differentiation import Grad
from ufl.core.multiindex import FixedIndex, MultiIndex

from dune.common.hashit import hashIt
from dune.source.builtin import get, hybridForEach, make_pair, make_index_sequence, make_shared
import dune.source.cplusplus as cplusplus
from dune.source.cplusplus import ConditionalExpression, Declaration, Using, Variable

from dune.source.cplusplus import AccessModifier, Declaration, Constructor, EnumClass, Include, InitializerList, Method, Struct, TypeAlias, UnformattedExpression, Variable
from dune.source.cplusplus import assign, construct, coordinate, dereference, lambda_, makeExpression, maxEdgeLength, minEdgeLength, return_

from .applyrestrictions import applyRestrictions

from dune.ufl.tensors import ExprTensor
from dune.ufl.tensors import apply as exprTensorApply
from dune.source.cplusplus import assign, construct, TypeAlias, Declaration, Variable,\
        UnformattedBlock, UnformattedExpression, Struct, return_,\
        SwitchStatement
from dune.source.cplusplus import Method as clsMethod
from dune.source.cplusplus import SourceWriter, ListWriter, StringWriter

from ufl import SpatialCoordinate,TestFunction,TrialFunction,Coefficient,\
        as_vector, as_matrix,dx,ds,grad,inner,zero,FacetNormal,dot
# from ufl.algorithms.analysis import extract_arguments_and_coefficients as coeff
from ufl.differentiation import Grad

def translateIndex(index):
    if isinstance(index, (tuple, MultiIndex)):
        return ''.join([translateIndex(i) for i in index])
    elif isinstance(index, (int, FixedIndex)):
        return '[ ' + str(index) + ' ]'
    else:
        raise Exception('Index type not supported: ' + repr(index))

class TooHighDerivative(Exception):
    '''codegen raises this when a too high derivative is encountered'''
    def __init__(self, error):
        Exception.__init__(self,error)

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

    # do nothing here (until complex conjugate is needed)
    def Conj(self,expr,x):
        return x

    # do nothing here (until complex conjugate is needed)
    def conj(self,expr,x):
        return x

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
                raise TooHighDerivative('CodeGenerator does not allow for third order derivatives, yet.')
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
        # set to False to disable tempVars
        # TODO: better dynamic switch
        # tempVars = False
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
            field = shape.ufl_function_space().field
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

def gridPartType(gf):
    try:
        gv = gf.space.grid._typeName
    except:
        gv = gf.grid._typeName
    gvType = re.split('::GridViewType$', gv)
    if len(gvType) == 2: # is a dune fem grid part
        return gvType[0]
    else:
        return 'Dune::FemPy::GridPart<'+gvType[0]+'>'


class ModelClass():
    def __init__(self, name, uflExpr, virtualize, dimRange=None, predefined=None):
        self.className = name
        self.targs = ['class GridPart']
        if dimRange is not None:
            self.bindable = True
            self.dimRange = dimRange
            self.bindableBase = 'Dune::Fem::BindableGridFunctionWithSpace<GridPart,Dune::Dim<'+str(self.dimRange)+'>>'
            self.bases = ['public '+self.bindableBase]
            self.includeFiles = ['dune/fem/function/localfunction/bindable.hh']
        else:
            self.bindable = False
            self.bases = []
            self.includeFiles = []
        self.includeFiles += ['dune/fem/common/intersectionside.hh']
        self.gridPartType = TypeAlias("GridPartType", "GridPart")
        self.ctor_args = []
        self.ctor_init = []
        self.skeleton = None
        self.init = None
        self.vars = None

        uflExpr = [e for e in uflExpr if e is not None]
        coefficients = set()
        for expr in uflExpr:
            try:
                coefficients |= set(expr.coefficients())
            except:
                _, cc = extract_arguments_and_coefficients(expr)
                coefficients |= set(cc)
        extracedAll = False
        while not extracedAll:
            extracedAll = True
            for c in coefficients:
                try:
                    predef = c.predefined
                except AttributeError:
                    continue
                for expr in predef.values():
                    _, cc = extract_arguments_and_coefficients(expr)
                    cc = set(cc)
                    if not cc.issubset(coefficients):
                       coefficients |= cc
                       extracedAll = False
                if not extracedAll:
                    break

        self.constantList = sorted((c for c in coefficients if c.is_cellwise_constant()), key=lambda c: c.count())
        self.coefficientList = sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())

        constants=[fieldVectorType(c,useScalar=True) for c in self.constantList]
        coefficients=(fieldVectorType(c) for c in self.coefficientList)
        constantNames=[getattr(c, 'name', None) for c in self.constantList]
        constantShapes=[getattr(c, 'ufl_shape', None) for c in self.constantList]
        coefficientNames=[getattr(c, 'name', None) for c in self.coefficientList]
        parameterNames=[getattr(c, 'parameter', None) for c in self.constantList]

        if not len(set(constantNames)) == len(constantNames):
            raise AttributeError("two constants have the same name which will lead to failure during code generation:"+','.join(c for c in constantNames))

        # this part modifies duplicate names in coefficient - slow version
        # can be improved
        coefficientNames_ = coefficientNames
        coefficientNames = []
        for c in coefficientNames_:
            if c is not None:
                while c in coefficientNames:
                    c = c+"A"
            coefficientNames.append(c)

        self._constants = [] if constants is None else list(constants)
        self._coefficients = [] if coefficients is None else list(coefficients)

        self._constantShapes = [None,] * len(self._constants) if constantShapes is None else list(constantShapes)
        self._constantNames  = [None,] * len(self._constants) if constantNames is None else list(constantNames)
        self._constantNames  = ['constant' + str(i) if n is None else n for i, n in enumerate(self._constantNames)]
        if len(self._constantNames) != len(self._constants):
            raise ValueError("Length of constantNames must match length of constants")
        invalidConstants = [n for n in self._constantNames if n is not None and re.match('^[a-zA-Z_][a-zA-Z0-9_]*$', n) is None]
        if invalidConstants:
            raise ValueError('Constant names are not valid C++ identifiers:' + ', '.join(invalidConstants) + '.')
        self._constantValues = [ None if not hasattr(c,"value") else c.value for c in self.constantList ]

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
        if self._coefficients:
            if virtualize:
                self.coefficientCppTypes = \
                    ['Dune::FemPy::VirtualizedGridFunction< ' +\
                     gridPartType(c) + ', ' + fieldVectorType(c) + ' >' \
                        if not c._typeName.startswith("Dune::Python::SimpleGridFunction") \
                        else c._typeName \
                    for c in self.coefficientList]
            else:
                self.coefficientCppTypes = [c._typeName for c in self.coefficientList]
        else:
            self.coefficientCppTypes = []

        # need to replace possible grid functions in values of predefined
        # import pdb; pdb.set_trace()
        self.predefined = {} if predefined is None else predefined
        for idx, coefficient in enumerate(self.coefficientList):
            try:
                self.predefined.update( coefficient.predefined )
            except AttributeError:
                pass
        for idx, coefficient in enumerate(self.coefficientList):
            for derivative in self.coefficient(idx, 'x', side=None):
                for k,v in self.predefined.items():
                    if coefficient == v:
                        self.predefined[k] = derivative
                coefficient = Grad(coefficient)

    @property
    def constantTypes(self):
        return ["Con" + n[0] + n[1:] for n in self._constantNames]

    @property
    def constantNames(self):
        return ["con" + n[0] + n[1:] for n in self._constantNames]

    @property
    def constantShortNames(self):
        return [n[0] + n[1:] for n in self._constantNames]

    @property
    def constantValues(self):
        return [str(0.) if v is None else\
                str(v) if not hasattr(v,"__len__") else\
                    str(list(v)).replace("[","{").replace("]","}")\
                for v in self._constantValues]

    @property
    def coefficientTypes(self):
        return ["Coeff" + n[0] + n[1:] for n in self._coefficientNames]

    @property
    def coefficientNames(self):
        return ["coeff" + n[0] + n[1:] for n in self._coefficientNames]

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
        predefined.update(self.predefined)

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


    def code(self, name=None, targs=None):
        # self.name = "Integrands"
        # self.targs = ['class GridPart']
        # self.gridPartType = TypeAlias("GridPartType", "GridPart")
        # self.ctor_args = []
        # self.ctor_init = []
        # self.skeleton (None is no intersecton terms are needed)
        # self.init (code to be added to init(Entity) method
        # self.vars (further class variables)
        # self.methods(code)

        if name is not None:
            self.name = name
        if targs is not None:
            self.targs = targs
        code = Struct(self.className,
                      targs=(self.targs + ['class ' + n for n in self.coefficientTypes]),
                      bases=self.bases)

        code.append(self.gridPartType)
        code.append(TypeAlias("GridView", "typename GridPartType::GridViewType"))
        if self.bindable:
            code.append(TypeAlias("FunctionSpaceType", "Dune::Fem::GridFunctionSpace<GridPartType,Dune::Dim<"+str(self.dimRange)+">>"))

        code.append(TypeAlias("EntityType", "typename GridPartType::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPartType::IntersectionType"))

        code.append(TypeAlias("GlobalCoordinateType", "typename EntityType::Geometry::GlobalCoordinate"))

        code.append(TypeAlias("Side","Dune::Fem::IntersectionSide"))

        for type, alias in zip(self._constants, self.constantTypes):
            code.append(TypeAlias(alias, type))
        constants = ["std::shared_ptr< " + c + " >" for c in self.constantTypes]
        if constants:
            code.append(TypeAlias("ConstantTupleType", "std::tuple< " + ", ".join(constants) + " >"))
            code.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantTupleType >::element_type", targs=["std::size_t i"]))
        else:
            code.append(TypeAlias("ConstantTupleType", "std::tuple<>"))

        if self._coefficients:
            coefficientSpaces = [('Dune::Fem::GridFunctionSpace< GridPartType, ' + c + ' >') for c in self._coefficients]
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " +", ".join(coefficientSpaces) + " >"))
            code.append(TypeAlias('CoefficientTupleType', 'std::tuple< ' + ', '.join(self.coefficientTypes) + ' >'))

            code.append(TypeAlias("CoefficientFunctionSpaceType", "std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >", targs=["std::size_t i"]))
            for s in ["RangeType", "JacobianRangeType"]:
                code.append(TypeAlias("Coefficient" + s, "typename CoefficientFunctionSpaceType< i >::" + s, targs=["std::size_t i"]))
        else:
            code.append(TypeAlias("CoefficientTupleType", "std::tuple<>"))

        code.append(TypeAlias('CoefficientType', 'std::tuple_element_t< i, CoefficientTupleType >', targs=['std::size_t i']))
        code.append(TypeAlias('ConstantType', 'typename std::tuple_element_t< i, ConstantTupleType >::element_type', targs=['std::size_t i']))

        # if self.skeleton is not None:
        #     code.append(EnumClass('Side', ['in = 0u', 'out = 1u'], 'std::size_t'))
        #     inside = '[ static_cast< std::size_t >( Side::in ) ]'
        # else:
        #     inside = ''

        if not self.bindable:
            if self.skeleton is None:
                entity_ = Variable('EntityType', 'entity_')
                insideEntity = entity_
            else:
                entity_ = Variable('std::array< EntityType, 2 >', 'entity_')
                insideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::in )')]
                outsideEntity = entity_[UnformattedExpression('std::size_t', 'static_cast< std::size_t >( Side::out )')]
            intersection_ = Variable('IntersectionType', 'intersection_')
        else:
            code.append(Using(self.bindableBase+'::entity'))

        constants_ = Variable("ConstantTupleType", "constants_")
        coefficientsTupleType = 'std::tuple< ' + ', '.join('Dune::Fem::ConstLocalFunction< ' + n + ' >' for n in self.coefficientTypes) + ' >'
        if self.skeleton is None:
            coefficients_ = Variable(coefficientsTupleType, 'coefficients_')
        else:
            coefficients_ = Variable('std::array< ' + coefficientsTupleType + ', 2 >', 'coefficients_')

        # generate code for constructor
        arg_param = Variable('const Dune::Fem::ParameterReader &', 'parameter')
        if self.bindable:
            args = [
                    Variable('const GridPartType &','gridPart'),
                    Variable('const std::string &','name'),
                    Variable('int','order')
                   ]
        else:
            args = []
        args += self.ctor_args + [Variable('const ' + t + ' &', n) for t, n in zip(self.coefficientTypes, self.coefficientNames)]
        if self.bindable:
            init = [self.bindableBase+'(gridPart,name,order)']
        else:
            init = []
        if self._coefficients:
            coeffInit = ['Dune::Fem::ConstLocalFunction< ' + n + ' >( ' + p + ' )' for n, p in zip(self.coefficientTypes, self.coefficientNames)]
            if self.skeleton is None:
                init += ["coefficients_( " + ", ".join(coeffInit) + " )"]
            else:
                init += ['coefficients_{{ ' + coefficientsTupleType + '( ' + ', '.join(coeffInit) + ' ), ' + coefficientsTupleType + '( ' + ', '.join(coeffInit) + ' ) }}']
        init = self.ctor_init + init
        args.append(Declaration(arg_param, initializer=UnformattedExpression('const ParameterReader &', 'Dune::Fem::Parameter::container()')))
        constructor = Constructor(args=args, init=init)
        for idx, (cppType, value) in enumerate(zip(self.constantTypes, self.constantValues)):
            constructor.append(assign(get(idx)(constants_),
                make_shared(cppType)(cppType+"(0)")))
        for idx, (name, cppType) in enumerate(zip(self._parameterNames, self.constantTypes)):
            if name is not None:
                constructor.append(assign(dereference(get(idx)(constants_)), UnformattedExpression('auto', arg_param.name + '.getValue< ' + cppType + ' >( "' + name + '" )', uses=[arg_param])))
        code.append(constructor)

        entity = Variable('const EntityType &', 'entity')
        intersection = Variable('const IntersectionType &', 'intersection')
        if self.bindable:
            initEntity = Method('void', 'bind', args=[entity])
            initEntity.append(self.bindableBase+'::bind(entity);')
            initIntersection = Method('void', 'bind', args=[intersection, Variable('Side', 'side')])
            initIntersection.append(self.bindableBase+'::bind(intersection,side);')
            code.append(initIntersection)
        else:
            initEntity = Method('bool', 'init', args=[entity])
            initEntity.append(assign(insideEntity, entity))
        if self.skeleton is None:
            for i, c in enumerate(self._coefficients):
                initEntity.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + ' ).bind( entity )', uses=[entity, coefficients_]))
        else:
            for i, c in enumerate(self._coefficients):
                initEntity.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::in ) ] ).bind( entity )', uses=[entity, coefficients_]))
        initEntity.append(self.init)
        if not self.bindable:
            initEntity.append(return_(True))
        code.append(initEntity)

        if not self.bindable:
            initIntersection = Method('bool', 'init', args=[intersection])
            initIntersection.append(assign(intersection_, intersection))
            if self.skeleton is None:
                initIntersection.append(return_('(intersection.boundary() && init( intersection.inside() ))'))
            else:
                initIntersection.append(assign(insideEntity, UnformattedExpression('EntityType', 'intersection.inside()')))
                for i, c in enumerate(self._coefficients):
                    # initIntersection.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::in ) ] ).bind( entity_[ static_cast< std::size_t >( Side::in ) ] )', uses=[coefficients_]))
                    initIntersection.append(UnformattedExpression('void', 'std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::in ) ] ).bind( intersection_, Side::in  )', uses=[coefficients_]))
                initIntersection.append('if( intersection.neighbor() )')
                initIntersection.append('{')
                initIntersection.append('  entity_[ static_cast< std::size_t >( Side::out ) ] = intersection.outside();')
                for i, c in enumerate(self._coefficients):
                    # initIntersection.append(UnformattedExpression('void', '  std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::out ) ] ).bind( entity_[ static_cast< std::size_t >( Side::out ) ] )', uses=[coefficients_]))
                    initIntersection.append(UnformattedExpression('void', '  std::get< ' + str(i) + ' >( ' + coefficients_.name + '[ static_cast< std::size_t >( Side::out ) ] ).bind( intersection_, Side::out )', uses=[coefficients_]))
                initIntersection.append('}')
                initIntersection.append(return_(True))
            code.append(initIntersection)

        ################################
        self.methods(code)
        ################################

        code.append(Method('const ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_))), const=True))
        code.append(Method('ConstantType< i > &', 'constant', targs=['std::size_t i'], code=return_(dereference(get('i')(constants_)))))

        for i, (t, n) in enumerate(zip(self.constantTypes, self.constantNames)):
            code.append(Method('const ' + t + ' &', n, code=return_(dereference(get(i)(constants_))), const=True))
            code.append(Method(t + ' &', n, code=return_(dereference(get(i)(constants_)))))

        if not self.bindable:
            code.append(Method('const EntityType &', 'entity', const=True, code=return_(insideEntity)))

        code.append(AccessModifier('private'))

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

        if not self.bindable:
            code.append(Declaration(entity_), Declaration(intersection_))
        code.append(Declaration(constants_), Declaration(coefficients_))

        if self.vars is not None:
            code += self.vars

        return code

def generateMethodBody(cppType, expr, returnResult, default, predefined):
    if expr is not None and not expr == [None]:
        try:
            dimR = expr.ufl_shape[0]
        except:
            if isinstance(expr,list) or isinstance(expr,tuple):
                expr = as_vector(expr)
            else:
                expr = as_vector([expr])
            dimR = expr.ufl_shape[0]

        _, coeff = extract_arguments_and_coefficients(expr)
        coeff = {c : c.toVectorCoefficient()[0] for c in coeff if len(c.ufl_shape) == 0 and not c.is_cellwise_constant()}
        expr = replace(expr, coeff)

        t = ExprTensor(expr.ufl_shape) # , exprTensorApply(lambda u: u, expr.ufl_shape, expr))
        expression = [expr[i] for i in t.keys()]
        u = extract_arguments_and_coefficients(expr)[0]
        if u != []:
            u = u[0]
            du = Grad(u)
            d2u = Grad(du)
            arg_u = Variable("const RangeType &", "u")
            arg_du = Variable("const JacobianRangeType &", "du")
            arg_d2u = Variable("const HessianRangeType &", "d2u")
            predefined.update( {u: arg_u, du: arg_du, d2u: arg_d2u} )
        code, results = generateCode(predefined, expression, tempVars=False)
        result = Variable(cppType, 'result')
        if cppType == 'double':
            code = code + [assign(result, results[0])]
        else:
            code = code + [assign(result[i], r) for i, r in zip(t.keys(), results)]
        if returnResult:
            code = [Declaration(result)] + code + [return_(result)]
    else:
        result = Variable(cppType, 'result')
        code = [assign(result, construct(cppType,default) )]
        if returnResult:
            code = [Declaration(result)] + code + [return_(result)]
    return code
def generateMethod(struct,expr, cppType, name,
        returnResult=True,
        defaultReturn='0',
        targs=None, args=None, static=False, const=False, volatile=False,
        inline=False,
        evalSwitch=True,
        predefined=None):
    if predefined is None:
        predefined = {}
    if not returnResult:
        args = args + [cppType + ' &result']
        returnType = 'void'
    else:
        returnType = cppType

    if isinstance(expr,dict):
        if evalSwitch:
            bndId = Variable('const int', 'bndId')
            code = SwitchStatement(bndId, default=return_(False))
            for id, e in expr.items():
                code.append(id,
                        [generateMethodBody('RangeType', e, False, defaultReturn,
                            predefined), return_(True)])
        else:
            code = UnformattedBlock()
        code = [code]
    else:
        code = generateMethodBody(cppType, expr, returnResult, defaultReturn, predefined)

    meth = clsMethod(returnType, name,
            code=code,
            args=args,
            targs=targs, static=static, const=const, volatile=volatile, inline=inline)
    struct.append(meth)

def uflSignature(form,*args):
    import ufl
    from ufl.algorithms.renumbering import renumber_indices
    sig = ''
    hashList  = [str(arg) for arg in args if not isinstance(arg,ufl.core.expr.Expr)]
    #   the following fails in the ufl algorithm in the rentrant corner problem:
    #   phi = atan_2(x[1], x[0]) + conditional(x[1] < 0, 2*math.pi, 0)
    #   create.function("ufl",gridview=grid,name="tmp",order=1,ufl=phi)
    # hashList += [str(renumber_indices(arg)) for arg in args if isinstance(arg,ufl.core.expr.Expr)]
    for arg in args:
        if isinstance(arg,ufl.core.expr.Expr):
            try:
                hashList += [str(renumber_indices(arg))]
            except ValueError:  # I don't think this should happen but it does :-<
                hashList += [str(arg)]

    if form is not None:
        hashList.append(form.signature())
    hashList = [h.split(" at ")[0] for h in hashList]
    return hashIt( sorted(hashList) )
