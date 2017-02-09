from __future__ import division, print_function, unicode_literals

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

from dune.source.builtin import get
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, dereference, nullptr, return_
from dune.source.cplusplus import ListWriter, SourceWriter
from dune.source.fem import declareFunctionSpace
from dune.generator import builder


class EllipticModel:
    def __init__(self, dimRange, signature):
        self.dimRange = dimRange
        self.coefficients = []
        self.init = None
        self.vars = None
        self.signature = signature
        self.field = "double"

        self.arg_r = Variable("RangeType &", "result")
        self.arg_dr = Variable("JacobianRangeType &", "result")

        self.arg_x = Variable("const Point &", "x")
        self.arg_u = Variable("const RangeType &", "u")
        self.arg_du = Variable("const JacobianRangeType &", "du")
        self.arg_d2u = Variable("const HessianRangeType &", "d2u")
        self.arg_ubar = Variable("const RangeType &", "ubar")
        self.arg_dubar = Variable("const JacobianRangeType &", "dubar")
        self.arg_d2ubar = Variable("const HessianRangeType &", "d2ubar")

        self.source = [assign(self.arg_r, construct("RangeType", 0))]
        self.linSource = [assign(self.arg_r, construct("RangeType", 0))]
        self.linNVSource = [assign(self.arg_r, construct("RangeType", 0))]
        self.diffusiveFlux = [assign(self.arg_dr, construct("JacobianRangeType", 0))]
        self.linDiffusiveFlux = [assign(self.arg_dr, construct("JacobianRangeType", 0))]
        self.fluxDivergence = [assign(self.arg_r, construct("RangeType", 0))]
        self.alpha = [assign(self.arg_r, construct("RangeType", 0))]
        self.linAlpha = [assign(self.arg_r, construct("RangeType", 0))]

        self.hasDirichletBoundary = False
        self.hasNeumanBoundary = False
        self.isDirichletIntersection = [return_(False)]
        self.dirichlet = [assign(self.arg_r, construct("RangeType", 0))]
        self.symmetric = False

    def pre(self, name):
        result = [TypeAlias("GridPartType", "GridPart"),
                  TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"),
                  TypeAlias("IntersectionType", "typename GridPart::IntersectionType")]

        result += declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), "GridPartType::dimensionworld", self.dimRange)
        result += [Declaration(Variable("const int", "dimLocal"), initializer="GridPartType::dimension", static=True)]

        constants = ["std::shared_ptr< Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " > >" for c in self.coefficients if c['constant']]
        result += [TypeAlias("ConstantsTupleType", "std::tuple< " + ", ".join(constants) + " >")]
        if constants:
            result += [TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantsTupleType >::element_type", targs=["std::size_t i"])]

        coefficients = ["Dune::Fem::FunctionSpace< DomainFieldType, " + SourceWriter.cpp_fields(c['field']) + ", dimDomain, " + str(c['dimRange']) + " >" for c in self.coefficients if not c['constant']]
        result += [Declaration(Variable("const std::size_t", "numCoefficients"), initializer=len(coefficients), static=True)]
        if coefficients:
            result += [TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " + ", ".join(coefficients) + " >")]
            result += [TypeAlias("Coefficient"+t, "typename std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >::" + t, targs=["std::size_t i"])
                  for t in ["RangeType", "JacobianRangeType", "HessianRangeType"]]

        result += [TypeAlias('CoefficientType', 'typename std::tuple_element< i, std::tuple< Coefficients... > >::type', targs=['std::size_t i'])]
        result += [TypeAlias('ConstantsType', 'typename std::tuple_element< i, ConstantsTupleType >::type::element_type', targs=['std::size_t i'])]

        result += [Constructor(code=UnformattedBlock("constructConstants( std::make_index_sequence< std::tuple_size<ConstantsTupleType>::value >() );"))]

        init = Method('bool', 'init', args=['const EntityType &entity'], const=True)
        init.append(UnformattedBlock("entity_ = &entity;", "initCoefficients( std::make_index_sequence< numCoefficients >() );"), self.init, return_(True))
        result += [init]

        result += [Method('const EntityType &', 'entity', code=return_(UnformattedExpression("const EntityType &", "*entity_")), const=True)]
        result += [Method('std::string', 'name', const=True, code=['return "' + name + '";'])]
        return result

    def main(self):
        result = [TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >")]
        result += [Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True)]

        result += [Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True),
                   Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True)]

        result += [Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True),
                   Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True)]

        result += [Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True)]

        result += [Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True),
                   Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True)]

        result += [Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)),
                   Method('bool', 'hasNeumanBoundary', const=True, code=return_(self.hasNeumanBoundary))]

        result += [Method('bool', 'isDirichletIntersection', args=['const IntersectionType &intersection', 'Dune::FieldVector< int, dimRange > &dirichletComponent'], code=self.isDirichletIntersection, const=True),
                   Method('void', 'dirichlet', targs=['class Point'], args=['int id', self.arg_x, self.arg_r], code=self.dirichlet, const=True)]
        return result

    def post(self):
        constants_ = Variable("ConstantsTupleType", "constants_")
        coefficients_ = Variable("std::tuple< Coefficients... >", "coefficients_")

        result = [Method("const ConstantsType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))), const=True),
                  Method("ConstantsType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))))]

        result += [Method("const CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)), const=True),
                   Method("CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)))]

        result += [AccessModifier("private")]

        initCoefficients = Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'], const=True)
        initCoefficients.append(UnformattedBlock("std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );"))

        constructConstants = Method('void', 'constructConstants', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        constructConstants.append(UnformattedBlock("std::ignore = std::make_tuple( (std::get< i >( constants_ ) = std::make_shared<ConstantsType< i >>(), i)... );"))

        result += [initCoefficients, constructConstants]

        result += [Declaration(Variable('const EntityType *', 'entity_'), nullptr, mutable=True),
                   Declaration(coefficients_, mutable=True), Declaration(constants_, mutable=True)]
        return result

    def code(self, name='Model', targs=[]):
        result = Struct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']))
        result.append(self.pre(name=name), self.main(), self.post())
        return result

    def write(self, sourceWriter, name='Model', targs=[]):
        sourceWriter.emit(self.code(name=name, targs=targs))

    def setCoef(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper'):
        sourceWriter.emit('')
        sourceWriter.typedef('std::tuple< ' + ', '.join(\
                [('Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< ' +\
                SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                str(coefficient['dimRange']) + ' > >') \
                for coefficient in self.coefficients if not coefficient["constant"]]) \
              + ' >', 'Coefficients')

        sourceWriter.openFunction('std::size_t renumberConstants', args=['pybind11::handle &obj'])
        sourceWriter.emit('std::string id = obj.str();')
        for coef in self.coefficients:
            number = str(coef['number'])
            name = coef['name']
            sourceWriter.emit('if (id == "' + name + '") return ' + number + ';')
        sourceWriter.emit('throw pybind11::value_error("coefficient \'" + id + "\' has not been registered");')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('void setConstant', targs=['std::size_t i'], args=[modelClass + ' &model', 'pybind11::list l'])
        sourceWriter.emit('model.template constant< i >() = l.template cast< typename ' + modelClass + '::ConstantsType<i> >();')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('auto defSetConstant', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.typedef('std::function< void( ' + modelClass + ' &model, pybind11::list ) >', 'Dispatch')
        sourceWriter.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setConstant< i > )... }};')
        sourceWriter.emit('')
        sourceWriter.emit('return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::list l ) {')
        sourceWriter.emit('    std::size_t k = renumberConstants(coeff);')
        sourceWriter.emit('    if( k >= dispatch.size() )')
        sourceWriter.emit('      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );' )
        sourceWriter.emit('    dispatch[ k ]( model.impl(), l );')
        sourceWriter.emit('    return k;')
        sourceWriter.emit('  };')
        sourceWriter.closeFunction()

        setCoefficient = Function('void', 'setCoefficient', targs=['std::size_t i'], args=[ modelClass + ' &model', 'pybind11::handle o'])
        setCoefficient.append('model.template coefficient< i >() = o.template cast< typename std::tuple_element< i, Coefficients >::type >().localFunction();')
        sourceWriter.emit(setCoefficient)

        defSetCoefficient = Function('auto', 'defSetCoefficient', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        defSetCoefficient.append(TypeAlias('Dispatch', 'std::function< void( ' + modelClass + ' &model, pybind11::handle ) >'),
                                 Declaration(Variable('std::array< Dispatch, sizeof...( i ) >', 'dispatch'), '{{ Dispatch( setCoefficient< i > )... }}'),
                                 '',
                                 'return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::handle o ) {',
                                 '    std::size_t k = renumberConstants(coeff);',
                                 'if( k >= dispatch.size() )',
                                 '      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );',
                                 '    dispatch[ k ]( model.impl(), o );',
                                 '    return k;'
                                 '  };')

        sourceWriter.emit(defSetCoefficient)

    def export(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper', constrArgs=(), constrKeepAlive=None):
        if self.coefficients:
            # sourceWriter.emit('cls.def( "setCoefficient", defSetCoefficient( std::make_index_sequence< std::tuple_size<Coefficients>::value >() ) );')
            sourceWriter.emit('cls.def( "setConstant", defSetConstant( std::make_index_sequence< std::tuple_size <typename '+ modelClass + '::ConstantsTupleType>::value >() ) );')
        sourceWriter.emit('')
        sourceWriter.emit('cls.def( "__init__", [] (' + wrapperClass + ' &instance, '+\
                ' '.join( i[1]+' '+i[0]+',' for i in constrArgs) +\
                'const pybind11::dict &coeff) {')
        sourceWriter.emit('  new (&instance) ' + wrapperClass + '('+\
                ', '.join(i[0] for i in constrArgs) +\
                ');')
        if self.coefficients:
            sourceWriter.emit('  const int size = std::tuple_size<Coefficients>::value;')
            sourceWriter.emit('  auto dispatch = defSetCoefficient( std::make_index_sequence<size>() );' )
            sourceWriter.emit('  std::vector<bool> coeffSet(size,false);')
            sourceWriter.emit('  for (auto item : coeff) {')
            sourceWriter.emit('    int k = dispatch(instance, item.first, item.second); ')
            sourceWriter.emit('    coeffSet[k] = true;')
            sourceWriter.emit('  }')
            sourceWriter.emit('  if ( !std::all_of(coeffSet.begin(),coeffSet.end(),[](bool v){return v;}) )')
            sourceWriter.emit('    throw pybind11::key_error("need to set all coefficients during construction");')
        sourceWriter.emit('  },')
        if constrKeepAlive:
            sourceWriter.emit(constrKeepAlive + ',')
        sourceWriter.emit(''.join('pybind11::arg("' + i[0] + '"), ' for i in constrArgs) +\
                'pybind11::arg("coefficients")=pybind11::dict() );')

    def codeCoefficient(self, code, coefficients, constants):
        """find coefficients/constants in code string and do replacements
        """
        if coefficients:
            number = 0
            numCoeffs = [coef['number'] for coef in self.coefficients if coef['constant'] == False]
            if numCoeffs:
                number = max(numCoeffs) + 1
            for key, val in coefficients.items():
                check = 0
                for coef in self.coefficients:
                    if key == coef['name']:
                        check = 1
                if check == 0:
                    self.coefficients.append({ \
                              'name' : key, \
                              'number' : number, \
                              'dimRange' : val.dimRange, \
                              'constant' : False, \
                              'field': "double" } )
                    number += 1
        if constants:
            number = 0
            numConsts = [coef['number'] for coef in self.coefficients if coef['constant'] == True]
            if numConsts:
                number = max(numCoeffs) + 1
            for key, val in constants.items():
                check = 0
                for coef in self.coefficients:
                    if key == coef['name']:
                        check = 1
                if check == 0:
                    self.coefficients.append({ \
                              'name' : key, \
                              'number' : number, \
                              'dimRange' : val, \
                              'constant' : True, \
                              'field': None } )
                    number += 1
        if '@const:' in code:
            codeCst = code.split('@const:')
            import itertools
            names = set(["".join(itertools.takewhile(str.isalpha, str(c))) for c in codeCst[1:]])
            number = 0
            numConsts = [coef['number'] for coef in self.coefficients if coef['constant'] == True]
            if numConsts:
                number = max(numConsts) + 1
            for name in names:
                check = 0
                for coef in self.coefficients:
                    if name == coef['name']:
                        check = 1
                if check == 0:
                    cname = '@const:' + name
                    afterName = code.split(cname)[1:]
                    if afterName[0][0] == '[':
                        beforeText = [an.split(']')[0].split('[')[1] for an in afterName]
                        dimRange = max( [int(bt) for bt in beforeText] ) + 1
                    else:
                        dimRange = 1
                    self.coefficients.append({
                          'name' : name, \
                          'number' : number, \
                          'dimRange' : dimRange, \
                          'constant' : True, \
                          'field': None } )
                    number += 1
        for coef in self.coefficients:
            num = str(coef['number'])
            gfname = '@gf:' + coef['name']
            constname = '@const:' + coef['name']
            jacname = '@jac:' + coef['name']
            if jacname in code:
                code = code.replace(jacname, 'dc' + num)
                if not 'CoefficientJacobianRangeType< ' + num + ' > dc' + num + ';' in code:
                    code = 'CoefficientJacobianRangeType< ' + num + ' > dc' + num + ';\n' \
                           + 'coefficient< ' + num + ' >().jacobian( x, dc' \
                           + num + ' );' + code
            if gfname in code:
                code = code.replace(gfname, 'c' + num)
                if not 'CoefficientRangeType< ' + num  + ' > c' + num + ';' in code:
                    code = 'CoefficientRangeType< ' + num  + ' > c' + num + ';\n' \
                           + 'coefficient< ' + num + ' >().evaluate( x, c' \
                           + num + ' );' + code
            elif constname in code:
                code = code.replace(constname, 'cc' + num)
                init = 'ConstantsRangeType< ' + num + ' > cc' + num + ' = constant< ' \
                           + num + ' >();'
                if not init in code:
                    code = init + code
        return code

    @staticmethod
    def codeDimRange(code):
        """find dimRange from code string
        """
        codeOut = ''
        lines = code.split("\n")
        if '@dimrange' in code or '@range' in code:
            for line in lines:
                if '@dimrange' in line or '@range' in line:
                    dimRange = int( line.split("=", 1)[1] )
                else:
                    codeOut += line + "\n"
            codeOut = codeOut[:-2]
        else:
            lhs = [line.split("=") for line in lines]
            values = [string[0] for string in lhs if "value" in string[0]]
            try:
                dimRange = max( [int(c.split("[")[1].split("]")[0]) for c in values] ) + 1
            except:
                dimRange = "n/a"
            codeOut = code
        return (codeOut, dimRange)

    def appendCode(self, key, code, **kwargs):
        function = getattr(self, key)
        coef = kwargs.pop("coefficients", {})
        const = kwargs.pop("constants", {})
        function.append(UnformattedBlock(self.codeCoefficient(code, coef, const)))
        setattr(self, key, function)



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

def generateCode(predefined, tensor, coefficients, tempVars=True):
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]



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


    arg = Variable('const RangeType &', 'u')
    darg = Variable('const JacobianRangeType &', 'du')
    d2arg = Variable('const HessianRangeType &', 'd2u')
    argbar = Variable('const RangeType &', 'ubar')
    dargbar = Variable('const JacobianRangeType &', 'dubar')
    d2argbar = Variable('const HessianRangeType &', 'd2ubar')
    model.source = generateCode({u: arg, du: darg, d2u: d2arg}, source, model.coefficients, tempVars)
    model.diffusiveFlux = generateCode({u: arg, du: darg}, diffusiveFlux, model.coefficients, tempVars)
    model.alpha = generateCode({u: arg}, boundarySource, model.coefficients, tempVars)
    model.linSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linSource, model.coefficients, tempVars)
    model.linNVSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linNVSource, model.coefficients, tempVars)
    model.linDiffusiveFlux = generateCode({u: arg, du: darg, ubar: argbar, dubar: dargbar}, linDiffusiveFlux, model.coefficients, tempVars)
    model.linAlpha = generateCode({u: arg, ubar: argbar}, linBoundarySource, model.coefficients, tempVars)
    model.fluxDivergence = generateCode({u: arg, du: darg, d2u: d2arg}, fluxDivergence, model.coefficients, tempVars)

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
            writer.emit(generateCode({}, ExprTensor((dimRange,), dirichlet[bndId]), model.coefficients, tempVars), indent=2, context=Method('void', 'dirichlet'))
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

    writer.emit(["#include <" + i + ">" for i in grid._includes])
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

    code = []

    nameSpace = NameSpace("ModelImpl_" + model.signature)
    nameSpace.append(model.code())
    code.append(nameSpace)

    code += [TypeAlias("GridPart", "typename Dune::FemPy::GridPart< " + grid._typeName + " >")]

    rangeTypes = ["Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " >" for c in model.coefficients if not c['constant']]
    coefficients = ["Dune::FemPy::VirtualizedLocalFunction< GridPart, " + r + " >" for r in rangeTypes]
    code += [TypeAlias("Model", nameSpace.name + "::Model< " + ", ".join(["GridPart"] + coefficients) + " >")]

    code += [TypeAlias("ModelWrapper", "DiffusionModelWrapper< Model >"),
             TypeAlias("ModelBase", "typename ModelWrapper::Base")]

    writer.emit(code)

    if model.coefficients:
        model.setCoef(writer)

    writer.openPythonModule(name)
    writer.emit('// export abstract base class')
    writer.emit('if( !pybind11::already_registered< ModelBase >() )')
    writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
    writer.emit('')
    writer.emit('// actual wrapper class for model derived from abstract base')
    writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
    writer.emit('')
    model.export(writer, 'Model', 'ModelWrapper')
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
