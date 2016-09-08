from __future__ import division, print_function

import ufl
import ufl.algorithms

import importlib
import os
import subprocess
import sys
import timeit
import types
from dune import comm
from dune.source import SourceWriter
from dune.source import BaseModel

# EllipticModel
# -------------

class EllipticModel(BaseModel):
    def __init__(self, dimRange, signature):
        BaseModel.__init__(self, dimRange, signature)
        self.source = "result = RangeType( 0 );"
        self.linSource = "result = RangeType( 0 );"
        self.linNVSource = "result = RangeType( 0 );"
        self.diffusiveFlux = "result = JacobianRangeType( 0 );"
        self.linDiffusiveFlux = "result = JacobianRangeType( 0 );"
        self.fluxDivergence = "result = RangeType( 0 );"
        self.alpha = "result = RangeType( 0 );"
        self.linAlpha = "result = RangeType( 0 );"
        self.hasDirichletBoundary = False
        self.hasNeumanBoundary = False
        self.isDirichletIntersection = "return false;"
        self.dirichlet = "result = RangeType( 0 );"
        self.f = "result = RangeType( 0 );"
        self.exact = "result = RangeType( 0 );"
        self.n = "result = RangeType( 0 );"
        self.jacobianExact = "result = JacobianRangeType( 0 );"

    def write(self, sourceWriter, name='Model', targs=[]):
        self.pre(sourceWriter, name='Model', targs=[])
        sourceWriter.typedef("Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >", "BoundaryIdProviderType")

        arg_x = 'const Point &x'
        arg_u = 'const RangeType &u'
        arg_du = 'const JacobianRangeType &du'
        arg_d2u = 'const HessianRangeType &d2u'
        arg_ubar = 'const RangeType &ubar'
        arg_dubar = 'const JacobianRangeType &dubar'
        arg_d2ubar = 'const HessianRangeType &d2ubar'
        arg_r = 'RangeType &result'
        arg_dr = 'JacobianRangeType &result'

        if True:
            sourceWriter.openConstMethod('void source', targs=['class Point'], args=[arg_x, arg_u, arg_du, arg_r])
            sourceWriter.emit(self.source)
            sourceWriter.closeConstMethod()

            sourceWriter.openConstMethod('void linSource', targs=['class Point'], args=[arg_ubar, arg_dubar, arg_x, arg_u, arg_du, arg_r])
            sourceWriter.emit(self.linSource)
            sourceWriter.closeConstMethod()
        else:
            sourceWriter.openConstMethod('void source', targs=['class Point'], args=[arg_x, arg_u, arg_du, arg_d2u, arg_r])
            sourceWriter.emit(self.source)
            sourceWriter.closeConstMethod()

            sourceWriter.openConstMethod('void linSource', targs=['class Point'], args=[arg_ubar, arg_dubar, arg_d2ubar, arg_x, arg_u, arg_du, arg_r])
            sourceWriter.emit(self.linSource)
            sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linNVSource', targs=['class Point'], args=[arg_ubar, arg_dubar, arg_d2ubar, arg_x, arg_d2u, arg_r])
        sourceWriter.emit(self.linNVSource)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void diffusiveFlux', targs=['class Point'], args=[arg_x, arg_u, arg_du, arg_dr])
        sourceWriter.emit(self.diffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linDiffusiveFlux', targs=['class Point'], args=[arg_ubar, arg_dubar, arg_x, arg_u, arg_du, arg_dr])
        sourceWriter.emit(self.linDiffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void fluxDivergence', targs=['class Point'], args=[arg_x, arg_u, arg_du, arg_d2u, arg_r])
        sourceWriter.emit(self.fluxDivergence)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void alpha', targs=['class Point'], args=[arg_x, arg_u, arg_r])
        sourceWriter.emit(self.alpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linAlpha', targs=['class Point'], args=[arg_ubar, arg_x, arg_u, arg_r])
        sourceWriter.emit(self.linAlpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasDirichletBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasDirichletBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasNeumanBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasNeumanBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool isDirichletIntersection', args=['const IntersectionType &intersection', 'Dune::FieldVector< int, dimRange > &dirichletComponent'])
        sourceWriter.emit(self.isDirichletIntersection)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void f', args=['const DomainType &x', arg_r])
        sourceWriter.emit(self.f)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void n', args=['const DomainType &x', arg_r])
        sourceWriter.emit(self.n)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void exact', args=['const DomainType &x', arg_r])
        sourceWriter.emit(self.exact)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void jacobianExact', args=['const DomainType &x', arg_dr])
        sourceWriter.emit('// used for possible computation of H^1 error')
        sourceWriter.emit(self.jacobianExact)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void dirichlet', targs=['class Point'], args=['int id', arg_x, arg_r])
        sourceWriter.emit(self.dirichlet)
        sourceWriter.closeConstMethod()

        self.post(sourceWriter, name='Model', targs=[])


# ExprTensor
# ----------

class ExprTensor:
    def __init__(self, shape, data=None):
        self.shape = shape
        if data is None:
            self.data = self._zero(shape)
        else:
            self.data = data

    def __repr__(self):
        return repr(self.data)

    def __add__(self, other):
        if not isinstance(other, ExprTensor):
            raise Exception('Cannot add ' + type(other) + ' to ' + type(self) + '.')
        if other.shape != self.shape:
            raise Exception('Cannot add tensors of different shape.')
        return ExprTensor(self.shape, self._add(self.shape, self.data, other.data))

    def __truediv__(self, other):
        if isinstance(other, ExprTensor):
            raise Exception('Cannot divide by tensors tensors.')
        return ExprTensor(self.shape, self._div(self.shape, self.data, other))

    def __mul__(self, other):
        if isinstance(other, ExprTensor):
            raise Exception('Cannot multiply tensors.')
        return ExprTensor(self.shape, self._mul(self.shape, self.data, other))

    def __getitem__(self, key):
        if isinstance(key, ufl.core.multiindex.MultiIndex):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(shape)) + '.')
            if not all(isinstance(idx, ufl.core.multiindex.FixedIndex) for idx in key):
                raise Exception('Expect only FixedIndex entries in MultiIndex.')
            data = self.data
            key = key.indices()
            while len(key) > 1:
                data = data[int(key[0])]
                key = key[1:]
            return data[int(key[0])]
        else:
            raise Exception('Expect MultiIndex to access component')

    def __setitem__(self, key, value):
        if isinstance(key, ufl.core.multiindex.MultiIndex):
            if len(key) != len(self.shape):
                raise Exception('Expect key of length ' + str(len(shape)) + '.')
            if not all(isinstance(idx, ufl.core.multiindex.FixedIndex) for idx in key):
                raise Exception('Expect only FixedIndex entries in MultiIndex.')
            data = self.data
            key = key.indices()
            while len(key) > 1:
                data = data[int(key[0])]
                key = key[1:]
            data[int(key[0])] = value
        else:
            raise Exception('Expect MultiIndex to access component')

    def keys(self):
        return [ufl.core.multiindex.MultiIndex(idx) for idx in self._keys(self.shape)]

    def as_ufl(self):
        return self._as_ufl(self.shape, self.data)

    def _add(self, shape, left, right):
        if len(shape) == 1:
            return [left[i] + right[i] for i in range(0, shape[0])]
        else:
            return [self._add(shape[1:], left[i], right[i]) for i in range(0, shape[0])]

    def _as_ufl(self, shape, data):
        if len(shape) == 1:
            return ufl.as_tensor(data)
        else:
            return ufl.as_tensor([self._as_ufl(shape[1:], data[i]) for i in range(0, shape[0])])

    def _keys(self, shape):
        if len(shape) == 1:
            return [(ufl.core.multiindex.FixedIndex(i),) for i in range(0, shape[0])]
        else:
            return [(ufl.core.multiindex.FixedIndex(i),) + t for i in range(0, shape[0]) for t in self._keys(shape[1:])]

    def _div(self, shape, tensor, value):
        if len(shape) == 1:
            return [tensor[i] / value for i in range(0, shape[0])]
        else:
            return [self._div(shape[1:], tensor[i], value) for i in range(0, shape[0])]

    def _mul(self, shape, tensor, value):
        if len(shape) == 1:
            return [tensor[i] * value for i in range(0, shape[0])]
        else:
            return [self._mul(shape[1:], tensor[i], value) for i in range(0, shape[0])]

    def _zero(self, shape):
        if len(shape) == 1:
            return [ufl.constantvalue.Zero() for i in range(0, shape[0])]
        else:
            return [self._zero(shape[1:]) for i in range(0, shape[0])]



# FluxExtracter
# -------------

class FluxExtracter(ufl.algorithms.transformer.Transformer):
    def __init__(self):
        ufl.algorithms.transformer.Transformer.__init__(self)

    def argument(self, expr):
        if expr.number() == 0:
            raise Exception('Test function should only occur in indexed expressions.')
        else:
            return expr

    grad = ufl.algorithms.transformer.Transformer.reuse_if_possible

    def indexed(self, expr):
        if len(expr.ufl_operands) != 2:
            raise Exception('Indexed expressions must have exactly two children.')
        operand = expr.ufl_operands[0]
        index = expr.ufl_operands[1]
        if self.isTestFunction(operand):
            tensor = ExprTensor(operand.ufl_shape)
            tensor[index] = ufl.constantvalue.IntValue(1)
            return {operand : tensor}
        else:
            return self.reuse_if_possible(expr, self.visit(operand), index)

    def division(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, ufl.core.expr.Expr):
            return {op : left[op] / right for op in left}
        else:
            raise Exception('Only the left child of a division may access the test function.')

    def product(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, ufl.core.expr.Expr):
            return {op : left[op] * right for op in left}
        elif isinstance(left, ufl.core.expr.Expr) and isinstance(right, dict):
            return {op : right[op] * left for op in right}
        else:
            raise Exception('Only one child of a product may access the test function.')

    def sum(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, dict):
            for op in right:
                left[op] = (left[op] + right[op]) if op in left else right[op]
            return left
        else:
            raise Exception('Either both summands must contain test function or none')

    def nonlinear(self, expr, *operands):
        for operand in operands:
            if isinstance(operand, dict):
                raise Exception('Test function cannot appear in nonlinear expressions.')
        return self.reuse_if_possible(expr, *operands)

    atan = nonlinear
    atan_2 = nonlinear
    cos = nonlinear
    sin = nonlinear
    power = nonlinear
    tan = nonlinear

    def terminal(self, expr):
        return expr

    def isTestFunction(self, expr):
        while isinstance(expr, ufl.differentiation.Grad):
            expr = expr.ufl_operands[0]
        return isinstance(expr, ufl.argument.Argument) and expr.number() == 0

    def variable(self, expr):
        return expr.expression()

# DerivativeExtracter
# -------------

class DerivativeExtracter(ufl.algorithms.transformer.Transformer):
    def __init__(self):
        ufl.algorithms.transformer.Transformer.__init__(self)

    def argument(self, expr):
        if expr.number() == 0:
            raise Exception('Test function should occur at all.')
        else:
            return expr

    grad = ufl.algorithms.transformer.Transformer.reuse_if_possible

    def indexed(self, expr):
        if len(expr.ufl_operands) != 2:
            raise Exception('Indexed expressions must have exactly two children.')
        operand = expr.ufl_operands[0]
        index = expr.ufl_operands[1]
        if self.isTrialFunction(operand):
            tensor = ExprTensor(operand.ufl_shape)
            tensor[index] = ufl.constantvalue.IntValue(1)
            return {operand : tensor}
        else:
            return self.reuse_if_possible(expr, self.visit(operand), index)

    def division(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, ufl.core.expr.Expr):
            return {op : left[op] / right for op in left}
        else:
            raise Exception('Only the left child of a division may access the test function.')

    def product(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, ufl.core.expr.Expr):
            return {op : left[op] * right for op in left}
        elif isinstance(left, ufl.core.expr.Expr) and isinstance(right, dict):
            return {op : right[op] * left for op in right}
        else:
            raise Exception('Only one child of a product may access the test function.')

    def sum(self, expr, left, right):
        if isinstance(left, ufl.core.expr.Expr) and isinstance(right, ufl.core.expr.Expr):
            return self.reuse_if_possible(expr, left, right)
        elif isinstance(left, dict) and isinstance(right, dict):
            for op in right:
                left[op] = (left[op] + right[op]) if op in left else right[op]
            return left
        else:
            raise Exception('Either both summands must contain test function or none')

    def nonlinear(self, expr, *operands):
        for operand in operands:
            if isinstance(operand, dict):
                raise Exception('Test function cannot appear in nonlinear expressions.')
        return self.reuse_if_possible(expr, *operands)

    atan = nonlinear
    atan_2 = nonlinear
    cos = nonlinear
    sin = nonlinear
    power = nonlinear
    tan = nonlinear

    def terminal(self, expr):
        return expr

    def isTrialFunction(self, expr):
        while isinstance(expr, ufl.differentiation.Grad):
            expr = expr.ufl_operands[0]
        return isinstance(expr, ufl.argument.Argument) and expr.number() == 1


# splitUFL2
# ------------
def splitUFL2(u,du,d2u,tree):
    tree0 = ExprTensor(u.ufl_shape)
    tree1 = ExprTensor(u.ufl_shape)
    tree2 = ExprTensor(u.ufl_shape)

    for index in tree.keys():
        q = DerivativeExtracter().visit(tree[index])
        for op in q:
            if op == u:
                tree0[index] = tree0[index] +\
                   sum(i[0]*i[1] for i in zip(q[op].as_ufl(),u))
            elif op == du:
                for r in range(du.ufl_shape[0]):
                    for d in range(du.ufl_shape[1]):
                        tree1[index] = tree1[index] +\
                            q[op].as_ufl()[r,d]*du[r,d]
            elif op == d2u:
                for r in range(d2u.ufl_shape[0]):
                    for d1 in range(d2u.ufl_shape[1]):
                        for d2 in range(d2u.ufl_shape[2]):
                            tree2[index] = tree2[index] +\
                                q[op].as_ufl()[r,d1,d2]*d2u[r,d1,d2]
            else:
                raise Exception('Invalid trial function derivative encountered in bulk integral: ' + op)
    return tree0,tree1,tree2

# splitUFLForm
# ------------

def splitUFLForm(form, linear):
    phi = form.arguments()[0]
    dphi = ufl.differentiation.Grad(phi)

    source = ExprTensor(phi.ufl_shape)
    diffusiveFlux = ExprTensor(dphi.ufl_shape)
    boundarySource = ExprTensor(phi.ufl_shape)

    form = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(form)))
    for integral in form.integrals():
        if integral.integral_type() == 'cell':
            fluxExprs = FluxExtracter().visit(integral.integrand())
            for op in fluxExprs:
                if op == phi:
                    source = source + fluxExprs[op]
                elif op == dphi:
                    diffusiveFlux = diffusiveFlux + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in bulk integral: ' + op)
        elif integral.integral_type() == 'exterior_facet':
            fluxExprs = FluxExtracter().visit(integral.integrand())
            for op in fluxExprs:
                if op == phi:
                    boundarySource = boundarySource + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in boundary integral: ' + op)
        else:
            raise NotImplementedError('Integrals of type ' + integral.integral_type() + ' are not supported.')

    if linear:
        u = form.arguments()[1]
        du = ufl.differentiation.Grad(u)
        d2u = ufl.differentiation.Grad(du)
        source0,source1,source2 = splitUFL2(u,du,d2u,source)
        source = source0 + source1 # + source2
        return source, source2, diffusiveFlux, boundarySource

    return source, diffusiveFlux, boundarySource

# CodeGenerator
# -------------

class CodeGenerator(ufl.algorithms.transformer.Transformer):
    def __init__(self, predefined, tempVars):
        ufl.algorithms.transformer.Transformer.__init__(self)
        self.using = set()
        self.exprs = predefined
        self.code = []
        self.tempVars = tempVars

    def argument(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            raise Exception('Unknown argument: ' + str(expr.number()))

    def atan(self, expr):
        self.using.add('using std::atan;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('atan( ' + self.visit(expr.ufl_operands[0]) + ' )')
        return self.exprs[expr]

    def atan_2(self, expr):
        self.using.add('using std::atan2;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('atan2( ' + self.visit(expr.ufl_operands[0]) + ', ' + self.visit(expr.ufl_operands[1]) + ' )')
        return self.exprs[expr]

    def coefficient(self, expr):
        if expr not in self.exprs:
            idx = str(expr.number)
            if expr.is_cellwise_constant():
                self.code.append('ConstantsRangeType< ' + idx + ' > cc' + idx + ' = constant< ' + idx + ' >();')
                self.exprs[expr] = 'cc' + idx
            else:
                self.code.append('CoefficientRangeType< ' + idx + ' > c' + idx + ';')
                self.code.append('coefficient< ' + idx + ' >().evaluate( x, c' + idx + ' );')
                self.exprs[expr] = 'c' + idx
        return self.exprs[expr]

    def cos(self, expr):
        self.using.add('using std::cos;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('cos( ' + self.visit(expr.ufl_operands[0]) + ' )')
        return self.exprs[expr]

    def division(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            left = self.visit(expr.ufl_operands[0])
            right = self.visit(expr.ufl_operands[1])
            self.exprs[expr] = self._makeTmp('(' + left + ' / ' + right + ')')
            return self.exprs[expr]

    def float_value(self, expr):
        val = str(expr.value())
        if "." not in val:
            val = val + "."
        if expr.value() < 0:
            return '(' + val + ')'
        else:
            return val

    def grad(self, expr):
        if expr not in self.exprs:
            operand = expr.ufl_operands[0]
            if isinstance(operand, ufl.coefficient.Coefficient):
                idx = str(operand.number)
                self.code.append('CoefficientJacobianRangeType< ' + idx + ' > dc' + idx + ';')
                self.code.append('coefficient< ' + idx + ' >().jacobian( x, dc' + idx + ' );')
                self.exprs[expr] = 'dc' + idx
            elif isinstance(operand, ufl.differentiation.Grad):
                operand = operand.ufl_operands[0]
                if isinstance(operand, ufl.coefficient.Coefficient):
                    idx = str(operand.number)
                    self.code.append('CoefficientHessianRangeType< ' + idx + ' > d2c' + idx + ';')
                    self.code.append('coefficient< ' + idx + ' >().hessian( x, d2c' + idx + ' );')
                    self.exprs[expr] = 'd2c' + idx
                else:
                    raise Exception('Elliptic model does not allow for second derivatives, yet.')
            elif isinstance(operand, ufl.argument.Argument):
                raise Exception('Unknown argument: ' + str(operand.number()))
            else:
                raise Exception('Cannot compute gradient of ' + repr(expr))
        return self.exprs[expr]

    def indexed(self, expr, operand, index):
        return operand + self.translateIndex(index)

    int_value = float_value

    def product(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            left = self.visit(expr.ufl_operands[0])
            right = self.visit(expr.ufl_operands[1])
            self.exprs[expr] = self._makeTmp('(' + left + ' * ' + right + ')')
            return self.exprs[expr]

    def power(self, expr):
        self.using.add('using std::pow;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('pow( ' + self.visit(expr.ufl_operands[0]) + ', ' + self.visit(expr.ufl_operands[1]) + ' )')
        return self.exprs[expr]

    def sin(self, expr):
        self.using.add('using std::sin;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('sin( ' + self.visit(expr.ufl_operands[0]) + ' )')
        return self.exprs[expr]

    def spatial_coordinate(self, expr):
        self.using.add('using Dune::Fem::coordinate;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('entity().geometry().global( coordinate( x ) )')
        return self.exprs[expr]

    # only implemented for 3D normal
    def facet_normal(self, expr):
        self.using.add('const DomainType w1 = -entity.geometry.jacobianTransposed(coordinate(x))[0];\n' \
                  '      const DomainType w2 = -entity.geometry.jacobianTransposed(coordinate(x))[1];\n' \
                  '      DomainType normal;\n' \
                  '      normal[0]=w1[1]*w2[2]-w1[2]*w2[1];\n' \
                  '      normal[1]=-(w1[0]*w2[2]-w1[2]*w2[0]);\n' \
                  '      normal[2]=w1[0]*w2[1]-w1[1]*w2[0];\n' \
                  '      normal/=2.*entity().geometry().volume();\n')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('normal')
        return self.exprs[expr]

    def sum(self, expr, left, right):
        return '(' + left + ' + ' + right + ')'

    def tan(self, expr):
        self.using.add('using std::tan;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('tan( ' + self.visit(expr.ufl_operands[0]) + ' )')
        return self.exprs[expr]

    def zero(self, expr):
        return '0'

    def _makeTmp(self, cexpr):
        if self.tempVars:
            var = 'tmp' + str(len(self.code))
            self.code.append('const auto ' + var + ' = ' + cexpr + ';')
            return var
        else:
            return cexpr

    def translateIndex(self, index):
        if isinstance(index, ufl.core.multiindex.MultiIndex):
            result = ''
            for component in index:
                result += self.translateIndex(component)
            return result
        elif isinstance(index, ufl.core.multiindex.FixedIndex):
            return '[ ' + str(index) + ' ]'
        else:
            raise Exception('Index type not supported: ' + repr(index))

# generateCode
# ------------

def generateCode(predefined, tensor, tempVars = True):
    generator = CodeGenerator(predefined, tempVars)
    results = []
    for index in tensor.keys():
        result = generator.visit(tensor[index])
        results.append('result' + generator.translateIndex(index) + ' = ' + result + ';')
    return list(generator.using) + generator.code + results



# compileUFL
# ----------

def compileUFL(equation, dirichlet = {}, exact = None, tempVars = True):
    form = equation.lhs - equation.rhs
    if not isinstance(form, ufl.Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    phi = form.arguments()[0]
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = ufl.differentiation.Grad(u)
    d2u = ufl.differentiation.Grad(du)
    ubar = ufl.Coefficient(u.ufl_function_space())
    dubar = ufl.differentiation.Grad(ubar)
    d2ubar = ufl.differentiation.Grad(dubar)

    field = u.ufl_function_space().ufl_element().field()

    # if exact solution is passed in subtract a(u,.) from the form
    if not exact == None:
        b = ufl.replace(form, {u: ufl.as_vector(exact)} )
        form = form - b

    dform = ufl.algorithms.apply_derivatives.apply_derivatives(ufl.derivative(ufl.action(form, ubar), ubar, u))

    source, diffusiveFlux, boundarySource = splitUFLForm( form, False )
    linSource, linNVSource, linDiffusiveFlux, linBoundarySource = splitUFLForm( dform, True )
    fluxDivergence, _, _ = splitUFLForm(ufl.inner(source.as_ufl() - ufl.div(diffusiveFlux.as_ufl()), phi) * ufl.dx(0),False)

    model = EllipticModel(dimRange, form.signature())

    model.field = field

    coefficients = set(form.coefficients())
    for bndId in dirichlet:
        for expr in dirichlet[bndId]:
            _, c = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
            coefficients |= set(c)

    idxConst = 0
    idxCoeff = 0
    for coefficient in coefficients:
        if coefficient.is_cellwise_constant():
            field = None  # must be improved for 'complex'
            idx = idxConst
            dimRange = 1 if coefficient.ufl_shape==() else coefficient.ufl_shape[0]
            idxConst += 1
        else:
            field = coefficient.ufl_function_space().ufl_element().field()
            dimRange = coefficient.ufl_shape[0]
            idx = idxCoeff
            dimRange = coefficient.ufl_shape[0]
            idxCoeff += 1
        setattr(coefficient,"number",idx)
        model.coefficients.append({ \
            'number' : coefficient.number, \
            'counter' : coefficient.count(), \
            'dimRange' : dimRange,\
            'constant' : coefficient.is_cellwise_constant(),\
            'field': field } )

    model.source = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, source, tempVars)
    model.diffusiveFlux = generateCode({ u : 'u', du : 'du' }, diffusiveFlux, tempVars)
    model.alpha = generateCode({ u : 'u' }, boundarySource, tempVars)
    model.linSource = generateCode({ u : 'u', du : 'du', d2u : 'd2u', ubar : 'ubar', dubar : 'dubar', d2ubar : 'd2ubar'}, linSource, tempVars)
    model.linNVSource = generateCode({ u : 'u', du : 'du', d2u : 'd2u', ubar : 'ubar', dubar : 'dubar', d2ubar : 'd2ubar'}, linNVSource, tempVars)
    model.linDiffusiveFlux = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linDiffusiveFlux, tempVars)
    model.linAlpha = generateCode({ u : 'u', ubar : 'ubar' }, linBoundarySource, tempVars)
    model.fluxDivergence = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, fluxDivergence, tempVars)

    if dirichlet:
        model.hasDirichletBoundary = True

        model.isDirichletIntersection = []
        model.isDirichletIntersection.append('const int bndId = BoundaryIdProviderType::boundaryId( intersection );')
        model.isDirichletIntersection.append('std::fill( dirichletComponent.begin(), dirichletComponent.end(), bndId );')
        model.isDirichletIntersection.append('switch( bndId )')
        model.isDirichletIntersection.append('{')
        for bndId in dirichlet:
            model.isDirichletIntersection.append('case ' + str(bndId) + ':')
        model.isDirichletIntersection.append('  return true;')
        model.isDirichletIntersection.append('default:')
        model.isDirichletIntersection.append('  return false;')
        model.isDirichletIntersection.append('}')

        model.dirichlet = []
        model.dirichlet.append('switch( id )')
        model.dirichlet.append('{')
        for bndId in dirichlet:
            if len(dirichlet[bndId]) != dimRange:
                raise Exception('Dirichtlet boundary condition has wrong dimension.')
            model.dirichlet.append('case ' + str(bndId) + ':')
            model.dirichlet.append('  {')
            model.dirichlet += ['    ' + line for line in generateCode({}, ExprTensor((dimRange,), dirichlet[bndId]), tempVars)]
            model.dirichlet.append('  }')
            model.dirichlet.append('  break;')
        model.dirichlet.append('default:')
        model.dirichlet.append('  result = RangeType( 0 );')
        model.dirichlet.append('}')

    return model



# importModel
# -----------

def importModel(grid, model, dirichlet = {}, exact = None, tempVars=True):
    start_time = timeit.default_timer()

    if isinstance(model, ufl.equation.Equation):
        model = compileUFL(model, dirichlet, exact, tempVars)
    compilePath = os.path.join(os.path.dirname(__file__), "../generated")

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + grid._typeHash

    if comm.rank == 0:
        print("Importing module for model with signature " + model.signature)
        if not os.path.isfile(os.path.join(compilePath, name + ".so")):
            writer = SourceWriter(compilePath + '/modelimpl.hh')
            writer.emit(grid._includes)
            writer.emit('')
            writer.emit('#include <dune/fem/misc/boundaryidprovider.hh>')
            writer.emit('#include <dune/fem/gridpart/leafgridpart.hh>')
            writer.emit('#include <dune/fem/gridpart/adaptiveleafgridpart.hh>')
            writer.emit('')
            writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
            writer.emit('#include <dune/corepy/pybind11/extensions.h>')
            writer.emit('')
            if model.coefficients:
                writer.emit('#include <dune/fempy/function/virtualizedgridfunction.hh>')
                writer.emit('')
            writer.emit('#include <dune/fem/schemes/diffusionmodel.hh>')

            modelNameSpace = 'ModelImpl_' + model.signature

            writer.openNameSpace(modelNameSpace)
            model.write(writer)
            writer.closeNameSpace(modelNameSpace)

            writer.typedef(grid._typeName, 'GridPart')

            if model.coefficients:
                writer.typedef(modelNameSpace + '::Model< GridPart' + ' '.join(\
                [(',Dune::FemPy::VirtualizedLocalFunction< GridPart,'+\
                    'Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' > >') \
                    for coefficient in model.coefficients if not coefficient["constant"]])\
                  + ' >', 'Model')
            else:
                writer.typedef(modelNameSpace + '::Model< GridPart >', 'Model')

            writer.typedef('DiffusionModelWrapper< Model >', 'ModelWrapper')
            writer.typedef('typename ModelWrapper::Base', 'ModelBase')

            if model.coefficients:
                model.set(writer)

            writer.openPythonModule(name)
            writer.emit('')
            writer.emit('')
            writer.emit('// export abstract base class')
            writer.emit('if( !pybind11::already_registered< ModelBase >() )')
            writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
            writer.emit('')
            writer.emit('// actual wrapper class for model derived from abstract base')
            writer.emit('pybind11::class_< ModelWrapper > model( module, "Model", pybind11::base< ModelBase >() );')
            writer.emit('model.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
            if model.coefficients:
                # index_sequence_for always gives a one element sequence (0)
                # writer.emit('model.def( "setCoefficient", defSetCoefficient( std::index_sequence_for< Coefficients >() ) );')
                writer.emit('model.def( "setCoefficient", defSetCoefficient( std::make_index_sequence< std::tuple_size<Coefficients>::value >() ) );')
                writer.emit('model.def( "setConstant", defSetConstant( std::make_index_sequence< std::tuple_size<typename Model::ConstantsTupleType>::value >() ) );')
            writer.emit('')
            writer.emit('module.def( "get", [] () { return new ModelWrapper(); } );')
            writer.closePythonModule(name)

            writer.close()

            # the new model is constructed in the file modelimpl.cc for which make targets exist:
            cmake = subprocess.Popen(["cmake", "--build", "../../..", "--target", "modelimpl"], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, "modelimpl.so"), os.path.join(compilePath, name + ".so"))
            print("Compilation took: " , timeit.default_timer()-start_time , "seconds")

        comm.barrier()
        return importlib.import_module("dune.generated." + name)
