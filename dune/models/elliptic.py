from __future__ import division, print_function

import ufl
import ufl.algorithms

import importlib
import os
import subprocess
import sys
import timeit
import types

from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices
from ufl.core.multiindex import FixedIndex, MultiIndex

from dune.ufl import codegen, GridCoefficient
from dune.ufl.tensors import ExprTensor
from dune.ufl.linear import splitMultiLinearExpr

from dune.source import Declaration, Method, TypeAlias, Variable
from dune.source.cplusplus import ListWriter, SourceWriter
from dune.source import BaseModel
from dune.generator import builder


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
        self.arg_x = 'const Point &x'
        self.arg_u = 'const RangeType &u'
        self.arg_du = 'const JacobianRangeType &du'
        self.arg_d2u = 'const HessianRangeType &d2u'
        self.arg_ubar = 'const RangeType &ubar'
        self.arg_dubar = 'const JacobianRangeType &dubar'
        self.arg_d2ubar = 'const HessianRangeType &d2ubar'
        self.arg_r = 'RangeType &result'
        self.arg_dr = 'JacobianRangeType &result'
        self.symmetric = 'false'

    def main(self, name='Model', targs=[]):
        hasDirichletBoundary = Method('bool', 'hasDirichletBoundary', const=True)
        hasDirichletBoundary.append('return ' + ('true' if self.hasDirichletBoundary else 'false') + ';')

        hasNeumanBoundary = Method('bool', 'hasNeumanBoundary', const=True)
        hasNeumanBoundary.append('return ' + ('true' if self.hasNeumanBoundary else 'false') + ';')

        isDirichletIntersection = Method('bool', 'isDirichletIntersection', args=['const IntersectionType &intersection', 'Dune::FieldVector< int, dimRange > &dirichletComponent'], code=self.isDirichletIntersection, const=True)

        dirichlet = Method('void', 'dirichlet', targs=['class Point'], args=['int id', self.arg_x, self.arg_r], code=self.dirichlet, const=True)

        result = []
        result.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
        result.append(Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True))

        result.append(Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True))
        result.append(Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True))

        result.append(Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True))
        result.append(Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True))

        result.append(Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True))

        result.append(Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True))
        result.append(Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True))

        result += [hasDirichletBoundary, hasNeumanBoundary, isDirichletIntersection, dirichlet]

        return result

    def write(self, sourceWriter, name='Model', targs=[]):
        self.pre(sourceWriter, name='Model', targs=[])
        sourceWriter.emit(self.main(name='Model', targs=[]))
        self.post(sourceWriter, name='Model', targs=[])

    def appendCode(self, key, code, **kwargs):
        function = getattr(self, key)
        newCode = '\n      '.join(function) + '\n' + code
        coef = kwargs.pop("coefficients", {})
        const = kwargs.pop("constants", {})
        newCode = self.codeCoefficient(newCode, coef, const)
        setattr(self, key, newCode)



# DerivativeExtracter
# -------------------

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
        if q==0: continue
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

    form = expand_indices(expand_derivatives(expand_compounds(form)))
    for integral in form.integrals():
        if integral.integral_type() == 'cell':
            fluxExprs = splitMultiLinearExpr(integral.integrand(), [phi])
            for op in fluxExprs:
                if op[0] == phi:
                    source = source + fluxExprs[op]
                elif op[0] == dphi:
                    diffusiveFlux = diffusiveFlux + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in bulk integral: ' + str(op[0]))
        elif integral.integral_type() == 'exterior_facet':
            fluxExprs = splitMultiLinearExpr(integral.integrand(), [phi])
            for op in fluxExprs:
                if op[0] == phi:
                    boundarySource = boundarySource + fluxExprs[op]
                else:
                    raise Exception('Invalid derivative encountered in boundary integral: ' + str(op[0]))
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



# generateCode
# ------------

def generateCode(predefined, tensor, coefficients, tempVars = True):
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars)
    return preamble + [('result' + codegen.translateIndex(i) + ' = ' + r + ';') for i, r in zip(keys, results)]



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

    try:
        field = u.ufl_function_space().ufl_element().field()
    except AttributeError:
        field = "double"

    # if exact solution is passed in subtract a(u,.) from the form
    if not exact == None:
        b = ufl.replace(form, {u: ufl.as_vector(exact)} )
        form = form - b

    dform = ufl.algorithms.apply_derivatives.apply_derivatives(ufl.derivative(ufl.action(form, ubar), ubar, u))

    source, diffusiveFlux, boundarySource = splitUFLForm( form, False )
    linSource, linNVSource, linDiffusiveFlux, linBoundarySource = splitUFLForm( dform, True )
    fluxDivergence, _, _ = splitUFLForm(ufl.inner(source.as_ufl() - ufl.div(diffusiveFlux.as_ufl()), phi) * ufl.dx(0),False)

    model = EllipticModel(dimRange, form.signature())

    model.hasNeumanBoundary = not boundarySource.is_zero()

    expandform = ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(equation.lhs)))
    if expandform == ufl.adjoint(expandform):
        model.symmetric = 'true'
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
            idxCoeff += 1
        try:
            name = coefficient.str()
        except:
            name = str(coefficient)
        model.coefficients.append({ \
            'name' : name, \
            'number' : idx, \
            'counter' : coefficient.count(), \
            'dimRange' : dimRange,\
            'constant' : coefficient.is_cellwise_constant(),\
            'field': field } )

    model.source = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, source, model.coefficients, tempVars)
    model.diffusiveFlux = generateCode({ u : 'u', du : 'du' }, diffusiveFlux, model.coefficients, tempVars)
    model.alpha = generateCode({ u : 'u' }, boundarySource, model.coefficients, tempVars)
    model.linSource = generateCode({ u : 'u', du : 'du', d2u : 'd2u', ubar : 'ubar', dubar : 'dubar', d2ubar : 'd2ubar'}, linSource, model.coefficients, tempVars)
    model.linNVSource = generateCode({ u : 'u', du : 'du', d2u : 'd2u', ubar : 'ubar', dubar : 'dubar', d2ubar : 'd2ubar'}, linNVSource, model.coefficients, tempVars)
    model.linDiffusiveFlux = generateCode({ u : 'u', du : 'du', ubar : 'ubar', dubar : 'dubar' }, linDiffusiveFlux, model.coefficients, tempVars)
    model.linAlpha = generateCode({ u : 'u', ubar : 'ubar' }, linBoundarySource, model.coefficients, tempVars)
    model.fluxDivergence = generateCode({ u : 'u', du : 'du', d2u : 'd2u' }, fluxDivergence, model.coefficients, tempVars)

    if dirichlet:
        model.hasDirichletBoundary = True

        writer = SourceWriter(ListWriter())
        writer.emit('const int bndId = BoundaryIdProviderType::boundaryId( intersection );')
        writer.emit('std::fill( dirichletComponent.begin(), dirichletComponent.end(), bndId );')
        writer.emit('switch( bndId )')
        writer.emit('{')
        for bndId in dirichlet:
            writer.emit('case ' + str(bndId) + ':')
        writer.emit('return true;', indent=1)
        writer.emit('default:')
        writer.emit('return false;', indent=1)
        writer.emit('}')
        model.isDirichletIntersection = writer.writer.lines

        writer = SourceWriter(ListWriter())
        writer.emit('switch( id )')
        writer.emit('{')
        for bndId in dirichlet:
            if len(dirichlet[bndId]) != dimRange:
                raise Exception('Dirichtlet boundary condition has wrong dimension.')
            writer.emit('case ' + str(bndId) + ':')
            writer.emit('{', indent=1)
            writer.emit(generateCode({}, ExprTensor((dimRange,), dirichlet[bndId]), model.coefficients, tempVars), indent=2)
            writer.emit('}', indent=1)
            writer.emit('break;', indent=1)
        writer.emit('default:')
        writer.emit('result = RangeType( 0 );', indent=1)
        writer.emit('}')
        model.dirichlet = writer.writer.lines

    return model



# generateModel
# -----------

def generateModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
    start_time = timeit.default_timer()

    if isinstance(model, ufl.equation.Equation):
        model = compileUFL(model, dirichlet, exact, tempVars)

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + grid._moduleName

    writer = SourceWriter()

    writer.emit("".join(["#include <" + i + ">\n" for i in grid._includes]))
    writer.emit('')
    writer.emit('#include <dune/fem/misc/boundaryidprovider.hh>')
    writer.emit('')
    writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
    writer.emit('#include <dune/corepy/pybind11/extensions.h>')
    writer.emit('')
    writer.emit('#include <dune/fempy/py/grid/gridpart.hh>')
    if model.coefficients:
        writer.emit('#include <dune/fempy/function/virtualizedgridfunction.hh>')
        writer.emit('')
    writer.emit('#include <dune/fem/schemes/diffusionmodel.hh>')

    modelNameSpace = 'ModelImpl_' + model.signature

    writer.openNameSpace(modelNameSpace)
    model.write(writer)
    writer.closeNameSpace(modelNameSpace)

    writer.typedef('typename Dune::FemPy::GridPart< ' + grid._typeName + ' >', 'GridPart')

    if model.coefficients:
        writer.typedef(modelNameSpace + '::Model< GridPart' + ' '.join(\
        [(', Dune::FemPy::VirtualizedLocalFunction< GridPart,'+\
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
        model.setCoef(writer)

    writer.openPythonModule(name)
    writer.emit('')
    writer.emit('')
    writer.emit('// export abstract base class')
    writer.emit('if( !pybind11::already_registered< ModelBase >() )')
    writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
    writer.emit('')
    writer.emit('// actual wrapper class for model derived from abstract base')
    writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
    writer.emit('')
    model.export(writer, 'Model', 'ModelWrapper')
    writer.emit('')
    writer.closePythonModule(name)

    if header != False:
        with open(header, 'w') as modelFile:
            modelFile.write(writer.writer.getvalue())
    return writer, name

def importModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
    if isinstance(model, str):
        with open(model, 'r') as modelFile:
            data = modelFile.read()
        name = data.split('PYBIND11_PLUGIN( ')[1].split(' )')[0]
        builder.load(name, data, "ellipticModel")
        return importlib.import_module("dune.generated." + name)
    writer, name = generateModel(grid, model, dirichlet, exact, tempVars, header)
    builder.load(name, writer.writer.getvalue(), "ellipticModel")
    writer.close()
    return importlib.import_module("dune.generated." + name)
