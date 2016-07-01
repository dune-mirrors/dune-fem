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
        elif isinstance(src, (list, tuple)):
            for srcline in src:
                self.emit(srcline)
        elif self._isstring(src):
            src = src.strip()
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

    def openStruct(self, name, targs=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + targs.strip() + ' >')
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
            self.emit('template< ' + targs.strip() + ' >')
        if args:
            self.emit(typedName + ' ( ' + args.strip() + ' ) const')
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
            self.emit('template< ' + targs.strip() + ' >')
        if args:
            self.emit(typedName + ' ( ' + args.strip() + ' )')
        else:
            self.emit(typedName + ' ()')
        self.emit('{')
        self.pushBlock('method', typedName)

    def closeMethod(self, typedName=None):
        self.popBlock('method', typedName)
        self.emit('}')

    def typedef(self, typeName, typeAlias, targs=None):
        if targs:
            self.emit('template< ' + targs.strip() + ' >')
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
    def __init__(self, dimRange):
        self.dimRange = dimRange
        self.init = None
        self.vars = None
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
        self.f = "result = RangeType( 0 );"
        self.g = "result = RangeType( 0 );"
        self.n = "result = RangeType( 0 );"
        self.jacobianExact = "result = JacobianRangeType( 0 );"

    def write(self, sourceWriter, name='Model', targs=None):
        if targs:
            targs = 'class GridPart, ' + targs.strip()
        else:
            targs = 'class GridPart'
        sourceWriter.openStruct(name, targs=targs)

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

        sourceWriter.openConstMethod('bool init', args='const EntityType &entity')
        sourceWriter.emit('entity_ = &entity;')
        sourceWriter.emit(self.init)
        sourceWriter.emit('return true;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('const EntityType &entity')
        sourceWriter.emit('return *entity_;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('std::string name')
        sourceWriter.emit('return "' + name + '";')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void source', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &result')
        sourceWriter.emit(self.source)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linSource', targs='class Point', args='const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &result')
        sourceWriter.emit(self.linSource)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void diffusiveFlux', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &result')
        sourceWriter.emit(self.diffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linDiffusiveFlux', targs='class Point', args='const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &result')
        sourceWriter.emit(self.linDiffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void fluxDivergence', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, const HessianRangeType &d2u, RangeType &result')
        sourceWriter.emit(self.fluxDivergence)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void alpha', targs='class Point', args='const Point &x, const RangeType &u, RangeType &result')
        sourceWriter.emit(self.alpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linAlpha', targs='class Point', args='const RangeType &ubar, const Point &x, const RangeType &u, RangeType &result')
        sourceWriter.emit(self.linAlpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasDirichletBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasDirichletBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasNeumanBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasNeumanBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool isDirichletIntersection', args='const IntersectionType &intersection, Dune::FieldVector< bool, dimRange > &dirichletComponent')
        sourceWriter.emit(self.isDirichletIntersection)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void f', args='const DomainType &x, RangeType &result')
        sourceWriter.emit(self.f)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void g', args='const DomainType &x, RangeType &result')
        sourceWriter.emit('// used both for dirichlet data and possible computation of L^2 error')
        sourceWriter.emit(self.g)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void n', args='const DomainType &x, RangeType &result')
        sourceWriter.emit(self.n)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void jacobianExact', args='const DomainType &x, JacobianRangeType &result')
        sourceWriter.emit('// used for possible computation of H^1 error')
        sourceWriter.emit(self.jacobianExact)
        sourceWriter.closeConstMethod()

        sourceWriter.section('private')
        sourceWriter.emit('mutable const EntityType *entity_ = nullptr;')
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

    def sin(self, expr, *operands):
        for operand in operands:
            if isinstance(operand, dict):
                raise Exception('Test function cannot appear in nonlinear expressions.')
        return self.reuse_if_possible(expr, *operands)

    cos = sin
    power = sin

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
    def __init__(self, predefined):
        ufl.algorithms.transformer.Transformer.__init__(self)
        self.using = set()
        self.exprs = predefined
        self.code = []

    def argument(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            raise Exception('Unknown argument: ' + str(expr.number()))

    def coefficient(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            raise Exception('Unknown coefficient: ' + repr(expr))

    def int_value(self, expr):
        if expr.value() < 0:
            return '(' + str(expr.value()) + ')'
        else:
            return str(expr.value())

    float_value = int_value

    def grad(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            operand = expr.ufl_operands[0]
            if isinstance(operand, ufl.differentiation.Grad):
                raise Exception('Elliptic model does not allow for second derivatives, yet.')
            elif isinstance(operand, ufl.argument.Argument):
                raise Exception('Unknown argument: ' + str(operand.number()))
            else:
                raise Exception('Cannot compute gradient of ' + repr(expr))

    def indexed(self, expr, operand, index):
        return operand + self.translateIndex(index)

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

    def sum(self, expr, left, right):
        return '(' + left + ' + ' + right + ')'

    def cos(self, expr):
        self.using.add('using std::cos;')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('cos( ' + self.visit(expr.ufl_operands[0]) + ' )')
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

    def zero(self, expr):
        return '0'

    def _makeTmp(self, cexpr):
        var = 'tmp' + str(len(self.code))
        self.code.append('const auto ' + var + ' = ' + cexpr + ';')
        return var

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

def generateCode(predefined, tensor):
    generator = CodeGenerator(predefined)
    results = []
    for index in tensor.keys():
        result = generator.visit(tensor[index])
        results.append('result' + generator.translateIndex(index) + ' = ' + result + ';')
    return list(generator.using) + generator.code + results



# compileUFL
# ----------

def compileUFL(equation):
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

    model = EllipticModel(dimRange)

    model.source = generateCode({ u : 'u', du : 'du' }, source)
    model.diffusiveFlux = generateCode({ u : 'u', du : 'du' }, diffusiveFlux)
    model.alpha = generateCode({ u : 'u' }, boundarySource)
    model.linSource = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linSource)
    model.linDiffusiveFlux = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linDiffusiveFlux)
    model.linAlpha = generateCode({ u : 'u', ubar : 'ubar' }, linBoundarySource)
    model.fluxDivergence = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, fluxDivergence)

    return model



# importModel
# -----------

def importModel(name, grid, model):
    compilePath = os.path.join(os.path.dirname(__file__), "../generated")

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'ellipticmodel_' + name + "_" + grid._typeHash

    if dune.femmpi.comm.rank == 0:
        writer = SourceWriter(compilePath + '/modelimpl.hh')
        writer.emit(grid._includes)
        for line in open(compilePath + '/modelimpl.hh.in', "rt"):
            if '#include "ModelTmp.hh"' in line:
                writer.openNameSpace('ModelTmp')
                model.write(writer)
                writer.closeNameSpace('ModelTmp')
            else:
              if '#MODELNAME' in line:
                  line = line.replace('#MODELNAME', name)
              if '#GRIDPARTCHOICE' in line:
                  line = line.replace('#GRIDPARTCHOICE', grid._typeName)
              if '#PYTEMPLATE' in line:
                  line = line.replace('#PYTEMPLATE', '')
              if '#DIMRANGE' in line:
                  line = line.replace('#DIMRANGE', str(model.dimRange))
              if '#PYSETCOEFFICIENT' in line:
                  line = ''
              elif '#PYRANGETYPE' in line:
                  line = line.replace('#PYRANGETYPE', '')
              writer.emit(line)
        writer.close()

        # the new model is constructed in the file modelimpl.cc for which make targets exist:
        cmake = subprocess.Popen(["cmake", "--build", "../../..", "--target", "modelimpl"], cwd=compilePath)
        cmake.wait()
        os.rename(os.path.join(compilePath, "modelimpl.so"), os.path.join(compilePath, name + ".so"))

        dune.femmpi.comm.barrier()
        return importlib.import_module("dune.generated." + name)
