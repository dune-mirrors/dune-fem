from __future__ import print_function

import ufl
import ufl.algorithms

# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")
        self.blocks = []
        self.begin = True

    def emit(self, src):
        if src is None:
            return
        elif isinstance(src, (list, tuple)):
            for srcline in src:
                self.emit(srcline)
        elif isinstance(src, str):
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



# EllipticModel
# -------------

class EllipticModel:
    def __init__(self, dimRange):
        self.dimRange = dimRange
        self.init = None
        self.vars = None
        self.source = "flux = RangeType( 0 );"
        self.linSource = "flux = JacobianRangeType( 0 );"
        self.diffusiveFlux = "flux = RangeType( 0 );"
        self.linDiffusiveFlux = "flux = JacobianRangeType( 0 );"
        self.fluxDivergence = "fluxDiv = RangeType( 0 );"
        self.alpha = "flux = RangeType( 0 );"
        self.linAlpha = "flux = RangeType( 0 );"
        self.hasDirichletBoundary = False
        self.hasNeumannBoundary = False
        self.isDirichletIntersection = "return false;"
        self.f = "value = RangeType( 0 );"
        self.g = "value = RangeType( 0 );"
        self.n = "value = RangeType( 0 );"
        self.jacobianExact = "value = JacobianRangeType( 0 );"

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
        sourceWriter.typedef("Dune::Fem::FunctionSpace< double, RangeFieldTaype, dimDomain, dimRange >", "FunctionSpaceType")
        sourceWriter.typedef("typename FunctionSpaceType::DomainType", "DomainType")
        sourceWriter.typedef("typename FunctionSpaceType::RangeType", "RangeType")
        sourceWriter.typedef("typename FunctionSpaceType::JacobianRangeType", "JacobianRangeType")
        sourceWriter.typedef("typename FunctionSpaceType::HessianRangeType", "HessianRangeType")

        sourceWriter.openConstMethod('bool init', args='const EntityType &entity')
        sourceWriter.emit('entity_ = &entity')
        sourceWriter.emit(self.init)
        sourceWriter.emit('return true;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('const EntityType &entity')
        sourceWriter.emit('return *entity_;')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('std::string name')
        sourceWriter.emit('return ' + name + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void source', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &flux')
        sourceWriter.emit(self.source)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linSource', targs='class Point', args='const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, RangeType &flux')
        sourceWriter.emit(self.linSource)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void diffusiveFlux', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &flux')
        sourceWriter.emit(self.diffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linDiffusiveFlux', targs='class Point', args='const RangeType &ubar, const JacobianRangeType &dubar, const Point &x, const RangeType &u, const JacobianRangeType &du, JacobianRangeType &flux')
        sourceWriter.emit(self.linDiffusiveFlux)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void fluxDivergence', targs='class Point', args='const Point &x, const RangeType &u, const JacobianRangeType &du, const HessianRangeType &d2u, RangeType &fluxDiv')
        sourceWriter.emit(self.fluxDivergence)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void alpha', targs='class Point', args='const Point &x, const RangeType &u, RangeType &flux')
        sourceWriter.emit(self.alpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void linAlpha', targs='class Point', args='const RangeType &ubar, const Point &x, const RangeType &u, RangeType &flux')
        sourceWriter.emit(self.linAlpha)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasDirichletBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasDirichletBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool hasNeumannBoundary')
        sourceWriter.emit('return ' + ('true' if self.hasNeumannBoundary else 'false') + ';')
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('bool isDirichletIntersection', args='const IntersectionType &intersection, Dune::FieldVector< bool, dimRange > &dirichletComponent')
        sourceWriter.emit(self.isDirichletIntersection)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void f', args='const DomainType &x, RangeType &value')
        sourceWriter.emit(self.f)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void g', args='const DomainType &x, RangeType &value')
        sourceWriter.emit('// used both for dirichlet data and possible computation of L^2 error')
        sourceWriter.emit(self.g)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void n', args='const DomainType &x, RangeType &value')
        sourceWriter.emit(self.n)
        sourceWriter.closeConstMethod()

        sourceWriter.openConstMethod('void jacobianExact', args='const DomainType &x, JacobianRangeType &jacobian')
        sourceWriter.emit('// used for possible computation of H^1 error')
        sourceWriter.emit(self.jacobianExact)
        sourceWriter.closeConstMethod()

        sourceWriter.section('private')
        sourceWriter.emit('mutable const EntityType *entity_ = nullptr;')
        sourceWriter.emit(self.vars)
        sourceWriter.closeStruct(name)



# FluxExtractor
# -------------

class FluxExtracter(ufl.algorithms.transformer.Transformer):
    def __init__(self):
        ufl.algorithms.transformer.Transformer.__init__(self)

    def argument(self, expr):
        if expr.number() == 0:
            return { 0 : expr }
        else:
            return { None : expr }

    def grad(self, expr, *children):
        if len(children) != 1:
            raise Exception('grad expressions must have exactly one child.')
        child = children[ 0 ]
        result = dict()
        for order in child:
            if order is None:
                result[None] = self.reuse_if_possible(expr, child[None])
            else:
                result[order+1] = self.reuse_if_possible(expr, child[order])
        return result

    def indexed(self, expr, *children):
        if len(children) != 2:
            raise Exception('indexed must have exactly two children.')
        if list(children[1].keys()) != [None]:
            raise Exception('second child of indexed cannot contain test function.')
        result = dict()
        for order in children[0]:
            result[order] = self.reuse_if_possible(expr, children[0][order], children[1][None])
        return result

    def product(self, expr, *children):
        testedChildren = self._testedChildren(*children)
        if not testedChildren:
            return { None : expr }
        elif len(testedChildren) == 1:
            result = dict()
            for order in testedChildren[0]:
                result[order] = self.reuse_if_possible(expr, *[child[order] if child == testedChildren[0] else child[None] for child in children])
            return result
        else:
            raise Exception('test function may not appear in nonlinear expression.')

    def sum(self, expr, *children):
        orders = set()
        for child in children:
            orders |= set(child.keys())
        result = dict()
        for order in orders:
            result[order] = self.reuse_if_possible(expr, *[child[order] if order in child else ufl.constantvalue.Zero() for child in children])
        return result

    def sin(self, expr, *children):
        if self._testedChildren(*children):
            raise Exception('test function may not appear in nonlinear expression.')
        return { None : expr }

    cos = sin

    def terminal(self, expr):
        return { None : expr }

    def _testedChildren(self, *children):
        result = []
        for child in children:
            if list(child.keys()) != [None]:
                result += [child]
        return result


# splitUFLForm
# ------------

def splitUFLForm(form):
    source = ufl.constantvalue.Zero()
    diffusiveFlux = ufl.constantvalue.Zero()

    form = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(form)))
    for integral in form.integrals():
        if integral.integral_type() == 'cell':
            splitIntegrands = FluxExtracter().visit(integral.integrand())
            if not (set(splitIntegrands.keys()) <= { 0, 1 }):
                raise Exception('Invalid derivatives encountered on test function: ' + repr(set(splitIntegrands.keys())))
            if 0 in splitIntegrands:
                source += splitIntegrands[0]
            if 1 in splitIntegrands:
                diffusiveFlux += splitIntegrands[1]
        else:
            raise NotImplementedError('Integrals of type ' + integral.integral_type() + ' are not supported.')

    return source, diffusiveFlux



# CodeGenerator
# -------------

class CodeGenerator(ufl.algorithms.transformer.Transformer):
    class Accessor:
        def __init__(self, index):
            self._index = index

        def index(self):
            return self._index

    def __init__(self, phi, predefined):
        ufl.algorithms.transformer.Transformer.__init__(self)
        self.phi = phi
        self.exprs = predefined
        self.code = []

    def argument(self, expr):
        if expr in self.exprs:
            return self.exprs[expr]
        elif expr == self.phi:
            raise Exception('Test function should only occur in indexed expressions.')
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
        if len(expr.ufl_operands) != 1:
            raise Exception('grad expressions must have exactly one child.')
        if expr in self.exprs:
            return self.exprs[expr]
        elif expr == self.phi:
            raise Exception('Test function should only occur in indexed expressions.')
        else:
            operand = expr.ufl_operands[0]
            if isinstance(operand, ufl.differentiation.Grad):
                raise Exception('Elliptic model does not allow for second derivatives, yet.')
            elif isinstance(operand, ufl.argument.Argument):
                raise Exception('Unknown argument: ' + str(operand.number()))
            else:
                raise Exception('Cannot compute gradient of ' + repr(expr))

    def indexed(self, expr):
        if len(expr.ufl_operands) != 2:
            raise Exception('indexed expressions must have exactly two children.')
        operand = expr.ufl_operands[0]
        index = expr.ufl_operands[1]
        if operand == self.phi:
            return CodeGenerator.Accessor(index)
        else:
            tensor = self.visit(operand)
            return tensor + self.translateIndex(index)

    def product(self, expr):
        if len(expr.ufl_operands) != 2:
            raise Exception('Product expressions must have exactly two children.')
        if expr in self.exprs:
            return self.exprs[expr]
        else:
            left = self.visit(expr.ufl_operands[0])
            right = self.visit(expr.ufl_operands[1])
            if not isinstance(left, str):
                left, right = right, left
            if not isinstance(left, str):
                raise Exception('Only one child of a product may access the test function.')
            if isinstance(right, CodeGenerator.Accessor):
                return { right.index() : left }
            if isinstance(right, dict):
                return { index : '(' + left + ' * ' + right[index] + ')' for index in right }
            else:
                self.exprs[expr] = self._makeTmp('(' + left + ' * ' + right + ')')
                return self.exprs[expr]

    def sum(self, expr):
        if len(expr.ufl_operands) != 2:
            raise Exception('Sum expressions must have exactly two children.')
        left = self.visit(expr.ufl_operands[0])
        right = self.visit(expr.ufl_operands[1])
        if isinstance(left, str) and isinstance(right, str):
            return '(' + left + ' + ' + right + ')'
        elif isinstance(left, dict) and isinstance(right, dict):
            for index in right:
                if index in left:
                    left[index] = '(' + left[index] + ' + ' + right[index] + ')'
                else:
                    left[index] = right[index]
            return left
        else:
            raise Exception('Either both summands must contain test function or none')

    def sin(self, expr):
        if len(expr.ufl_operands) != 1:
            raise Exception('sin expressions must have exactly one child.')
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('std::sin( ' + self.visit(expr.ufl_operands[0]) + ' )')
        return self.exprs[expr]

    def spatial_coordinate(self, expr):
        if expr not in self.exprs:
            self.exprs[expr] = self._makeTmp('entity().geometry().global( x )')
        return self.exprs[expr]

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
            return '[' + str(index) + ']'
        else:
            raise Exception('Index type not supported: ' + repr(index))




# generateCode
# ------------

def generateCode(phi, predefined, expr):
    generator = CodeGenerator(phi, predefined)
    result = generator.visit(expr)
    code = generator.code
    if not isinstance(result, dict):
        raise Exception('Expression did not contain test function.')
    for index in result:
        code.append('flux' + generator.translateIndex(index) + ' = ' + result[index] + ';')
    return code



# compileUFL
# ----------

def compileUFL(equation, dimRange):
    form = equation.lhs - equation.rhs
    if not isinstance(form, ufl.Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires from with at least two arguments.")

    source, diffusiveFlux = splitUFLForm( form )

    phi = form.arguments()[0]
    dphi = ufl.differentiation.Grad(phi)
    u = form.arguments()[1]
    du = ufl.differentiation.Grad(u)
    ubar = ufl.Coefficient(u.ufl_function_space())
    dubar = ufl.differentiation.Grad(ubar)
    dform = ufl.algorithms.apply_derivatives.apply_derivatives(ufl.derivative(ufl.action(form, ubar), ubar, u))

    linSource, linDiffusiveFlux = splitUFLForm( dform )

    model = EllipticModel(dimRange)

    model.source = generateCode(phi, { u : 'u', du : 'du' }, source)
    model.diffusiveFlux = generateCode(dphi, { u : 'u', du : 'du' }, diffusiveFlux)
    model.linSource = generateCode(phi, { u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linSource)
    model.linDiffusiveFlux = generateCode(dphi, { u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linDiffusiveFlux)

    return model
