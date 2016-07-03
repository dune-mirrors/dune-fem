from __future__ import print_function

import ufl
import ufl.algorithms

import importlib
import os
import subprocess
import sys
import types
import dune.femmpi


# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")
        self.blocks = []
        self.begin = True

    def close(self):
        if self.blocks:
            raise Exception("Open blocks left in source file.")
        self.file.close()

    def emit(self, src):
        if src is None:
            return
        elif isinstance(src, (list, set, tuple)):
            for srcline in src:
                self.emit(srcline)
        elif self._isstring(src):
            src = src.rstrip()
            if src:
                print('  ' * len(self.blocks) + src, file=self.file)
            else:
                print(file=self.file)
            self.begin = False
        else:
            raise Exception("Unable to print " + repr(src) + ".")

    def pushBlock(self, block, name):
        self.blocks.append((block, name))
        self.begin = True

    def popBlock(self, block, name=None):
        (realBlock, realName) = self.blocks.pop()
        if realBlock != block or (name and name != realName):
            self.blocks.append((realBlock, realName))
            raise Exception("Trying to close " + realBlock + " " + realName + " as " + block + " " + name + ".");

    def openNameSpace(self, name):
        self.emit(None if self.begin else '')
        self.emit('namespace ' + name)
        self.emit('{')
        self.emit('')
        self.pushBlock('namespace', name)

    def closeNameSpace(self, name=None):
        self.emit(None if self.begin else '')
        self.popBlock('namespace', name)
        self.emit('} // namespace ' + name)

    def openClass(self, name, targs=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('class ' + name)
        self.emit('{')
        self.pushBlock('class', name)

    def closeClass(self, name=None):
        self.popBlock('class', name)
        self.emit('};')

    def openStruct(self, name, targs=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('struct ' + name)
        self.emit('{')
        self.pushBlock('struct', name)

    def closeStruct(self, name=None):
        self.popBlock('struct', name)
        self.emit('};')

    def section(self, section):
        self.emit(None if self.begin else '')
        (block, name) = self.blocks.pop()
        if block != 'class' and block != 'struct':
            self.blocks.push((block, name))
            raise Exception("Trying to declare " + section.strip() + " section outside class or struct.")
        self.emit(section.strip() + ':')
        self.blocks.append((block, name))
        self.begin = True

    def openConstMethod(self, typedName, targs=None, args=None):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' ) const')
        else:
            self.emit(typedName + ' () const')
        self.emit('{')
        self.pushBlock('const method', typedName)

    def closeConstMethod(self, typedName=None):
        self.popBlock('const method', typedName)
        self.emit('}')

    def openMethod(self, typedName, targs=None, args=None):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' )')
        else:
            self.emit(typedName + ' ()')
        self.emit('{')
        self.pushBlock('method', typedName)

    def closeMethod(self, typedName=None):
        self.popBlock('method', typedName)
        self.emit('}')

    def openFunction(self, typedName, targs=None, args=None):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' )')
        else:
            self.emit(typedName + ' ()')
        self.emit('{')
        self.pushBlock('function', typedName)

    def closeFunction(self, typedName=None):
        self.popBlock('function', typedName)
        self.emit('}')

    def openPythonModule(self, moduleName):
        self.emit(None if self.begin else '')
        self.emit('PYBIND11_PLUGIN( ' + moduleName.strip() + ' )')
        self.emit('{')
        self.pushBlock('pybind11 module', moduleName)
        self.emit('pybind11::module module( "' + moduleName.strip() + '" );')

    def closePythonModule(self, moduleName=None):
        self.emit('')
        self.emit('return module.ptr();')
        self.popBlock('pybind11 module', moduleName)
        self.emit('}')

    def typedef(self, typeName, typeAlias, targs=None):
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
            self.emit('using ' + typeAlias + ' = ' + typeName + ';')
        else:
            self.emit('typedef ' + typeName + ' ' + typeAlias + ';')

    def _isstring(self, obj):
        if sys.version_info.major == 2:
            return isinstance(obj, basestring)
        else:
            return isinstance(obj, str)



# EllipticModel
# -------------

class EllipticModel:
    def __init__(self, dimRange, signature):
        self.dimRange = dimRange
        self.coefficients = []
        self.init = None
        self.vars = None
        self.signature = signature
        self.source = "result = RangeType( 0 );"
        self.linSource = "result = JacobianRangeType( 0 );"
        self.diffusiveFlux = "result = RangeType( 0 );"
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
        sourceWriter.openStruct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']))

        sourceWriter.typedef("GridPart", "GridPartType")
        sourceWriter.typedef("double", "RangeFieldType")

        sourceWriter.emit("static const int dimRange = " + str(self.dimRange) + ";")
        sourceWriter.emit("static const int dimDomain = GridPartType::dimensionworld;")
        sourceWriter.emit("static const int dimLocal = GridPartType::dimension;")

        sourceWriter.typedef("typename GridPart::template Codim< 0 >::EntityType", "EntityType")
        sourceWriter.typedef("typename GridPart::IntersectionType", "IntersectionType")
        sourceWriter.typedef("Dune::Fem::FunctionSpace< double, RangeFieldType, dimDomain, dimRange >", "FunctionSpaceType")
        sourceWriter.typedef("typename FunctionSpaceType::DomainType", "DomainType")
        sourceWriter.typedef("typename FunctionSpaceType::RangeType", "RangeType")
        sourceWriter.typedef("typename FunctionSpaceType::JacobianRangeType", "JacobianRangeType")
        sourceWriter.typedef("typename FunctionSpaceType::HessianRangeType", "HessianRangeType")

        sourceWriter.typedef("Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >", "BoundaryIdProviderType")

        sourceWriter.section('private')
        if self.coefficients:
            sourceWriter.typedef('std::tuple< ' + ', '.join([('Dune::Fem::FunctionSpace< double, RangeFieldType, dimDomain, ' + str(coefficient['dimRange']) + ' >') for coefficient in self.coefficients]) + ' >', 'CoefficientFunctionSpaceTupleType')
            sourceWriter.section('public')
            sourceWriter.emit('static const std::size_t numCoefficients = std::tuple_size< CoefficientFunctionSpaceTupleType >::value;')
            sourceWriter.typedef('typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type', 'CoefficientFunctionSpaceType', targs=['std::size_t i'] )
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::RangeType', 'CoefficientRangeType', targs=['std::size_t i'])
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::JacobianRangeType', 'CoefficientJacobianRangeType', targs=['std::size_t i'])
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::HessianRangeType', 'CoefficientHessianRangeType', targs=['std::size_t i'])
        else:
            sourceWriter.emit('static const std::size_t numCoefficients = 0u;')
            sourceWriter.section('public')

        sourceWriter.emit('')
        sourceWriter.typedef('typename std::tuple_element< i, std::tuple< Coefficients... > >::type', 'CoefficientType', targs=['std::size_t i'])

        arg_x = 'const Point &x'
        arg_u = 'const RangeType &u'
        arg_du = 'const JacobianRangeType &du'
        arg_d2u = 'const HessianRangeType &d2u'
        arg_ubar = 'const RangeType &ubar'
        arg_dubar = 'const JacobianRangeType &dubar'
        arg_r = 'RangeType &result'
        arg_dr = 'JacobianRangeType &result'

        sourceWriter.openConstMethod('bool init', args=['const EntityType &entity'])
        sourceWriter.emit('entity_ = &entity;')
        sourceWriter.emit('initCoefficients( std::make_index_sequence< numCoefficients >() );')
        sourceWriter.emit(self.init)
        sourceWriter.emit('return true;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('const EntityType &entity')
        sourceWriter.emit('return *entity_;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('std::string name')
        sourceWriter.emit('return "' + name + '";')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void source', targs=['class Point'], args=[arg_x, arg_u, arg_du, arg_r])
        sourceWriter.emit(self.source)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linSource', targs=['class Point'], args=[arg_ubar, arg_dubar, arg_x, arg_u, arg_du, arg_r])
        sourceWriter.emit(self.linSource)
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

        sourceWriter.openConstMethod('const CoefficientType< i > &coefficient', targs=['std::size_t i'])
        sourceWriter.emit('return std::get< i >( coefficients_ );')
        sourceWriter.closeConstMethod()
        sourceWriter.openMethod('CoefficientType< i > &coefficient', targs=['std::size_t i'])
        sourceWriter.emit('return std::get< i >( coefficients_ );')
        sourceWriter.closeMethod()

        sourceWriter.section('private')
        sourceWriter.openConstMethod('void initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.emit('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );')
        sourceWriter.closeConstMethod()
        sourceWriter.emit('')
        sourceWriter.emit('mutable const EntityType *entity_ = nullptr;')
        sourceWriter.emit('mutable std::tuple< Coefficients... > coefficients_;')
        sourceWriter.emit(self.vars)
        sourceWriter.closeStruct(name)



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
            raise Exception('indexed expressions must have exactly two children.')
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



# splitUFLForm
# ------------

def splitUFLForm(form):
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
            idx = str(expr.count())
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
        if expr.value() < 0:
            return '(' + str(expr.value()) + ')'
        else:
            return str(expr.value())

    def grad(self, expr):
        if expr not in self.exprs:
            operand = expr.ufl_operands[0]
            if isinstance(operand, ufl.coefficient.Coefficient):
                idx = str(operand.count())
                self.code.append('CoefficientJacobianRangeType< ' + idx + ' > dc' + idx + ';')
                self.code.append('coefficient< ' + idx + ' >().jacobian( x, dc' + idx + ' );')
                self.exprs[expr] = 'dc' + idx
            elif isinstance(operand, ufl.differentiation.Grad):
                operand = operand.ufl_operands[0]
                if isinstance(operand, ufl.coefficient.Coefficient):
                    idx = str(operand.count())
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

def compileUFL(equation, dirichlet = {}, tempVars = True):
    form = equation.lhs - equation.rhs
    if not isinstance(form, ufl.Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires from with at least two arguments.")

    phi = form.arguments()[0]
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = ufl.differentiation.Grad(u)
    d2u = ufl.differentiation.Grad(du)
    ubar = ufl.Coefficient(u.ufl_function_space())
    dubar = ufl.differentiation.Grad(ubar)

    dform = ufl.algorithms.apply_derivatives.apply_derivatives(ufl.derivative(ufl.action(form, ubar), ubar, u))

    source, diffusiveFlux, boundarySource = splitUFLForm( form )
    linSource, linDiffusiveFlux, linBoundarySource = splitUFLForm( dform )
    fluxDivergence, _, _ = splitUFLForm(ufl.inner(- ufl.div(diffusiveFlux.as_ufl()), phi) * ufl.dx(0))

    model = EllipticModel(dimRange, form.signature())

    for coefficient in form.coefficients():
        model.coefficients.append({'dimRange' : coefficient.ufl_shape[0]})

    model.source = generateCode({ u : 'u', du : 'du' }, source, tempVars)
    model.diffusiveFlux = generateCode({ u : 'u', du : 'du' }, diffusiveFlux, tempVars)
    model.alpha = generateCode({ u : 'u' }, boundarySource, tempVars)
    model.linSource = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linSource, tempVars)
    model.linDiffusiveFlux = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linDiffusiveFlux, tempVars)
    model.linAlpha = generateCode({ u : 'u', ubar : 'ubar' }, linBoundarySource, tempVars)
    model.fluxDivergence = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, fluxDivergence, tempVars)

    if dirichlet:
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

def importModel(grid, model, dirichlet = {}):
    if isinstance(model, ufl.equation.Equation):
        model = compileUFL(model, dirichlet)
    compilePath = os.path.join(os.path.dirname(__file__), "../generated")

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + grid._typeHash

    if dune.femmpi.comm.rank == 0:
        print("Importing module for model with signature " + model.signature)
        if not os.path.isfile(os.path.join(compilePath, name + ".so")):
            writer = SourceWriter(compilePath + '/modelimpl.hh')
            writer.emit(grid._includes)
            writer.emit('')
            writer.emit('#include <dune/fem/gridpart/leafgridpart.hh>')
            writer.emit('#include <dune/fem/gridpart/adaptiveleafgridpart.hh>')
            writer.emit('')
            writer.emit('#include <dune/fempy/pybind11/pybind11.h>')
            writer.emit('#include <dune/fempy/pybind11/extensions.h>')
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
                writer.typedef(modelNameSpace + '::Model< GridPart, ' + ', '.join([('Dune::FemPy::VirtualizedLocalFunction< GridPart, Dune::FieldVector< double, ' + str(coefficient['dimRange']) + ' > >') for coefficient in model.coefficients])  + ' >', 'Model')
            else:
                writer.typedef(modelNameSpace + '::Model< GridPart >', 'Model')

            writer.typedef('DiffusionModelWrapper< Model >', 'ModelWrapper')
            writer.typedef('typename ModelWrapper::Base', 'ModelBase')

            if model.coefficients:
                writer.emit('')
                writer.typedef('std::tuple< ' + ', '.join([('Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, ' + str(coefficient['dimRange']) + ' > >') for coefficient in model.coefficients]) + ' >', 'Coefficients')

                writer.openFunction('void setCoefficient', targs=['std::size_t i'], args=['Model &model', 'pybind11::object o'])
                writer.emit('model.template coefficient< i >() = o.template cast< typename std::tuple_element< i, Coefficients >::type >().localFunction();')
                writer.closeFunction()

                writer.openFunction('auto defSetCoefficient', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
                writer.typedef('std::function< void( Model &model, pybind11::object ) >', 'Dispatch')
                writer.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setCoefficient< i > )... }};')
                writer.emit('')
                writer.emit('return [ dispatch ] ( ModelWrapper &model, std::size_t k, pybind11::object o ) {')
                writer.emit('    if( k >= dispatch.size() )')
                writer.emit('      throw std::range_error( "No such coefficient." );')
                writer.emit('    dispatch[ k ]( model.impl(), o );')
                writer.emit('  };')
                writer.closeFunction()

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
                writer.emit('model.def( "setCoefficient", defSetCoefficient( std::index_sequence_for< Coefficients >() ) );')
            writer.emit('')
            writer.emit('module.def( "get", [] () { return new ModelWrapper(); } );')
            writer.closePythonModule(name)

            writer.close()

            # the new model is constructed in the file modelimpl.cc for which make targets exist:
            cmake = subprocess.Popen(["cmake", "--build", "../../..", "--target", "modelimpl"], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, "modelimpl.so"), os.path.join(compilePath, name + ".so"))

        dune.femmpi.comm.barrier()
        return importlib.import_module("dune.generated." + name)
