from __future__ import division, print_function

from ufl.corealg.map_dag import map_expr_dags
from ufl.corealg.multifunction import MultiFunction
from ufl.argument import Argument
from ufl.coefficient import Coefficient
from ufl.differentiation import Grad
from ufl.core.multiindex import FixedIndex, MultiIndex

import dune.source.cplusplus as cplusplus
from dune.source.cplusplus import ConditionalExpression, Declaration, Using, Variable

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

    def cos(self, expr, x):
        self.using.add(Using(cplusplus.cos))
        return self._makeTmp(cplusplus.cos(x))

    def division(self, expr, x, y):
        return self._makeTmp(x / y)

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

    def lt(self, expr, left, right):
        return left < right

    def max_value(self, expr, left, right):
        self.using.add(Using(cplusplus.max_))
        return self._makeTmp(cplusplus.max_(left, right))

    def min_value(self, expr, left, right):
        self.using.add(Using(cplusplus.min_))
        return self._makeTmp(cplusplus.min_(left, right))

    def multi_index(self, expr):
        return expr

    int_value = float_value

    def negative_restricted(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        operand = expr.ufl_operands[0]
        raise Exception('Cannot compute restriction of ' + str(operand) + ', yet')

    def positive_restricted(self, expr):
        try:
            return self._makeTmp(self.predefined[expr], True)
        except KeyError:
            pass

        operand = expr.ufl_operands[0]
        raise Exception('Cannot compute restriction of ' + str(operand) + ', yet')

    def product(self, expr, x, y):
        return self._makeTmp(x * y)

    def power(self, expr, x, y):
        self.using.add(Using(cplusplus.pow_))
        return self._makeTmp(cplusplus.pow_(x, y))

    def sin(self, expr, x):
        self.using.add(Using(cplusplus.sin))
        return self._makeTmp(cplusplus.sin(x))

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

    def zero(self, expr):
        return cplusplus.makeExpression(0)

    def _getNumber(self, expr):
        try:
            name = expr.str()
        except:
            name = str(expr)
        e = [ee for ee in self.coefficients if ee["name"] == name]
        if len(e) > 1:
            raise KeyError('two coefficients provided with same name')
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
    generator = CodeGenerator(predefined, coefficients, tempVars)
    results = map_expr_dags(generator, expressions)
    return list(generator.using) + generator.code, results
