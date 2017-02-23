from __future__ import division, print_function, unicode_literals

import importlib
import os
import subprocess
import sys
import timeit
import types

from ufl import Coefficient, Form, SpatialCoordinate
from ufl import action, adjoint, derivative, div, dx, inner
from ufl.algorithms import expand_compounds, expand_derivatives, expand_indices
from ufl.algorithms.analysis import extract_arguments_and_coefficients
from ufl.algorithms.apply_derivatives import apply_derivatives
from ufl.differentiation import Grad
from ufl.equation import Equation
from ufl.core.multiindex import FixedIndex, MultiIndex

from dune.ufl import DirichletBC, GridCoefficient
from dune.ufl import codegen
from dune.ufl.tensors import ExprTensor
from dune.ufl.linear import splitMultiLinearExpr

from dune.source.builtin import get, hybridForEach, make_index_sequence, make_shared
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, SwitchStatement, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, dereference, lambda_, nullptr, return_
from dune.source.cplusplus import ListWriter, SourceWriter
from dune.source.fem import declareFunctionSpace
from dune.generator import builder
from dune.common.hashit import hashIt


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

        self.arg_i = Variable("const IntersectionType &", "intersection")
        self.arg_bndId = Variable("int", "bndId")

        self.source = [assign(self.arg_r, construct("RangeType", 0))]
        self.linSource = [assign(self.arg_r, construct("RangeType", 0))]
        # self.linNVSource = [assign(self.arg_r, construct("RangeType", 0))]
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

    def code(self, name='Model', targs=[]):
        code = Struct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']))

        code.append(TypeAlias("GridPartType", "GridPart"))
        code.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))

        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), "GridPartType::dimensionworld", self.dimRange))
        code.append(Declaration(Variable("const int", "dimLocal"), initializer="GridPartType::dimension", static=True))

        constants = ["std::shared_ptr< Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " > >" for c in self.coefficients if c['constant']]
        code.append(TypeAlias("ConstantsTupleType", "std::tuple< " + ", ".join(constants) + " >"))
        if constants:
            code.append(TypeAlias("ConstantsRangeType", "typename std::tuple_element_t< i, ConstantsTupleType >::element_type", targs=["std::size_t i"]))

        coefficients = ["Dune::Fem::FunctionSpace< DomainFieldType, " + SourceWriter.cpp_fields(c['field']) + ", dimDomain, " + str(c['dimRange']) + " >" for c in self.coefficients if not c['constant']]
        code.append(Declaration(Variable("const std::size_t", "numCoefficients"), initializer=len(coefficients), static=True))
        if coefficients:
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " + ", ".join(coefficients) + " >"))
            code.append([TypeAlias("Coefficient"+t, "typename std::tuple_element_t< i, CoefficientFunctionSpaceTupleType >::" + t, targs=["std::size_t i"])
                         for t in ["RangeType", "JacobianRangeType", "HessianRangeType"]])

        code.append(TypeAlias('CoefficientType', 'typename std::tuple_element< i, std::tuple< Coefficients... > >::type', targs=['std::size_t i']))
        code.append(TypeAlias('ConstantsType', 'typename std::tuple_element< i, ConstantsTupleType >::type::element_type', targs=['std::size_t i']))

        entity_ = Variable('const EntityType *', 'entity_')
        constants_ = Variable("ConstantsTupleType", "constants_")
        coefficients_ = Variable("std::tuple< Coefficients... >", "coefficients_")

        code.append(Constructor(code=hybridForEach(make_index_sequence("std::tuple_size< ConstantsTupleType >::value")(),
                                                   lambda_(args=["auto i"], capture=[Variable('auto', 'this')], code=assign(get("i")(constants_), make_shared("ConstantsType< i >")())))))

        init = Method('bool', 'init', args=['const EntityType &entity'], const=True)
        init.append(UnformattedBlock("entity_ = &entity;", "initCoefficients( std::make_index_sequence< numCoefficients >() );"), self.init, return_(True))
        code.append(init)

        code.append(Method('const EntityType &', 'entity', code=return_(dereference(entity_)), const=True))
        code.append(Method('std::string', 'name', const=True, code=return_('"' + name + '"')))

        code.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
        code.append(Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True))

        code.append(Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True))
        code.append(Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True))

        code.append(Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True))
        code.append(Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True))

        code.append(Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True))

        code.append(Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True))
        code.append(Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True))

        code.append(Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)))
        code.append(Method('bool', 'hasNeumanBoundary', const=True, code=return_(self.hasNeumanBoundary)))

        code.append(Method('bool', 'isDirichletIntersection', args=[self.arg_i, 'Dune::FieldVector< int, dimRange > &dirichletComponent'], code=self.isDirichletIntersection, const=True))
        code.append(Method('void', 'dirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_r], code=self.dirichlet, const=True))

        code.append(Method("const ConstantsType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))), const=True))
        code.append(Method("ConstantsType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_)))))

        code.append(Method("const CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)), const=True))
        code.append(Method("CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_))))

        code.append(AccessModifier("private"))

        initCoefficients = Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'], const=True)
        initCoefficients.append(UnformattedBlock("std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );"))

        code.append(initCoefficients)

        code.append(Declaration(entity_, nullptr, mutable=True))
        code.append(Declaration(coefficients_, mutable=True), Declaration(constants_, mutable=True))
        return code

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



# splitUFLForm
# ------------

def splitUFLForm(form):
    phi = form.arguments()[0]
    dphi = Grad(phi)

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

    return source, diffusiveFlux, boundarySource



# splitUFL2
# ---------

def splitUFL2(u,du,d2u,tree):
    tree0 = ExprTensor(u.ufl_shape)
    tree1 = ExprTensor(u.ufl_shape)
    tree2 = ExprTensor(u.ufl_shape)

    for index in tree.keys():
        q = splitMultiLinearExpr(tree[index], [u])
        if q==0: continue
        for op in q:
            if not isinstance(op, tuple) or (len(op) != 1):
                raise Exception('Missing trial function in bulk integral')
            if op[0] == u:
                tree0[index] += sum(i[0]*i[1] for i in zip(q[op].as_ufl(),u))
            elif op[0] == du:
                for r in range(du.ufl_shape[0]):
                    for d in range(du.ufl_shape[1]):
                        tree1[index] += q[op].as_ufl()[r,d]*du[r,d]
            elif op[0] == d2u:
                for r in range(d2u.ufl_shape[0]):
                    for d1 in range(d2u.ufl_shape[1]):
                        for d2 in range(d2u.ufl_shape[2]):
                            tree2[index] += q[op].as_ufl()[r,d1,d2]*d2u[r,d1,d2]
            else:
                raise Exception('Invalid trial function derivative encountered in bulk integral: ' + str(op[0]))
    return tree0, tree1, tree2



# generateCode
# ------------

def generateCode(predefined, tensor, coefficients, tempVars=True):
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]


#def compileUFL(equation, dirichlet = {}, exact = None, tempVars = True):
def compileUFL(equation, *args, **kwargs):
    form = equation.lhs - equation.rhs
    if not isinstance(form, Form):
        raise Exception("ufl.Form expected.")
    if len(form.arguments()) < 2:
        raise Exception("Elliptic model requires form with at least two arguments.")

    phi = form.arguments()[0]
    dimRange = phi.ufl_shape[0]

    u = form.arguments()[1]
    du = Grad(u)
    d2u = Grad(du)
    ubar = Coefficient(u.ufl_function_space())
    dubar = Grad(ubar)
    d2ubar = Grad(dubar)

    x = SpatialCoordinate(form.ufl_cell())

    try:
        field = u.ufl_function_space().ufl_element().field()
    except AttributeError:
        field = "double"

    # if exact solution is passed in subtract a(u,.) from the form
    if "exact" in kwargs:
        b = ufl.replace(form, {u: ufl.as_vector(kwargs["exact"])} )
        form = form - b

    dform = apply_derivatives(derivative(action(form, ubar), ubar, u))

    source, diffusiveFlux, boundarySource = splitUFLForm(form)
    linSource, linDiffusiveFlux, linBoundarySource = splitUFLForm(dform)
    fluxDivergence, _, _ = splitUFLForm(inner(source.as_ufl() - div(diffusiveFlux.as_ufl()), phi) * dx(0))

    # split linNVSource off linSource
    # linSources = splitUFL2(u, du, d2u, linSource)
    # linNVSource = linSources[2]
    # linSource = linSources[0] + linSources[1]

    model = EllipticModel(dimRange, form.signature())

    model.hasNeumanBoundary = not boundarySource.is_zero()

    expandform = expand_indices(expand_derivatives(expand_compounds(equation.lhs)))
    if expandform == adjoint(expandform):
        model.symmetric = 'true'
    model.field = field

    dirichletBCs = [arg for arg in args if isinstance(arg, DirichletBC)]
    if "dirichlet" in kwargs:
        dirichletBCs += [DirichletBC(u.ufl_space(), ufl.as_vector(value), bndId) for bndId, value in kwargs["dirichlet"].items()]

    coefficients = set(form.coefficients())
    for bc in dirichletBCs:
        _, c = extract_arguments_and_coefficients(bc.value)
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

    tempVars = kwargs.get("tempVars", True)

    predefined = {u: model.arg_u, du: model.arg_du}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.source = generateCode(predefined, source, model.coefficients, tempVars)
    model.diffusiveFlux = generateCode(predefined, diffusiveFlux, model.coefficients, tempVars=tempVars)
    predefined.update({ubar: model.arg_ubar, dubar: model.arg_dubar})
    model.linSource = generateCode(predefined, linSource, model.coefficients, tempVars=tempVars)
    model.linDiffusiveFlux = generateCode(predefined, linDiffusiveFlux, model.coefficients, tempVars=tempVars)

    # model.linNVSource = generateCode({u: arg, du: darg, d2u: d2arg, ubar: argbar, dubar: dargbar, d2ubar: d2argbar}, linNVSource, model.coefficients, tempVars)

    predefined = {u: model.arg_u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.alpha = generateCode(predefined, boundarySource, model.coefficients, tempVars)
    predefined.update({ubar: model.arg_ubar})
    model.linAlpha = generateCode(predefined, linBoundarySource, model.coefficients, tempVars)

    predefined = {u: model.arg_u, du: model.arg_du, d2u: model.arg_d2u}
    predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
    model.fluxDivergence = generateCode(predefined, fluxDivergence, model.coefficients, tempVars=tempVars)

    if dirichletBCs:
        model.hasDirichletBoundary = True

        bySubDomain = dict()
        for bc in dirichletBCs:
            if bc.subDomain in bySubDomain:
                raise Exception('Multiply defined Dirichlet boundary for subdomain ' + str(bc.subDomain))

            if not isinstance(bc.functionSpace, (FunctionSpace, FiniteElementBase)):
                raise Exception('Function space must either be a ufl.FunctionSpace or a ufl.FiniteElement')
            if isinstance(bc.functionSpace, FunctionSpace) and (bc.functionSpace != u.ufl_function_space()):
                raise Exception('Cannot handle boundary conditions on subspaces, yet')
            if isinstance(bc.functionSpace, FiniteElementBase) and (bc.functionSpace != u.ufl_element()):
                raise Exception('Cannot handle boundary conditions on subspaces, yet')

            value = ExprTensor(u.ufl_shape)
            for key in value.keys():
                value[key] = Indexed(bc.value, MultiIndex(tuple(FixedIndex(k) for k in key)))
            bySubDomain[bc.subDomain] = value

        bndId = Variable('const int', 'bndId')
        getBndId = UnformattedExpression('int', 'Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >::boundaryId( ' + model.arg_i.name + ' )')

        switch = SwitchStatement(bndId, default=return_(False))
        for i in bySubDomain:
            switch.append(i, return_(True))
        model.isDirichletIntersection = [Declaration(bndId, initializer=UnformattedExpression('int', 'BoundaryIdProviderType::boundaryId( ' + model.arg_i.name + ' )')),
                                         UnformattedBlock('std::fill( dirichletComponent.begin(), dirichletComponent.end(), ' + bndId.name + ' );'),
                                         switch
                                        ]

        switch = SwitchStatement(model.arg_bndId, default=assign(model.arg_r, construct("RangeType", 0)))
        predefined = {}
        predefined[x] = UnformattedExpression('auto', 'entity().geometry().global( Dune::Fem::coordinate( ' + model.arg_x.name + ' ) )')
        for i, value in bySubDomain.items():
            switch.append(i, generateCode(predefined, value, model.coefficients, tempVars=tempVars))
        model.dirichlet = [switch]

    return model


#def generateModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
def generateModel(grid, model, *args, **kwargs):
    start_time = timeit.default_timer()

    if isinstance(model, Equation):
        model = compileUFL(model, *args, **kwargs)

    # if not isinstance(grid, types.ModuleType):
    #     grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + hashIt(grid._typeName)

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

    if "header" in kwargs:
        with open(kwargs["header"], 'w') as modelFile:
            modelFile.write(writer.writer.getvalue())
    return writer, name


#def importModel(grid, model, dirichlet = {}, exact = None, tempVars = True, header = False):
def importModel(grid, model, *args, **kwargs):
    if isinstance(model, str):
        with open(model, 'r') as modelFile:
            data = modelFile.read()
        name = data.split('PYBIND11_PLUGIN( ')[1].split(' )')[0]
        builder.load(name, data, "ellipticModel")
        return importlib.import_module("dune.generated." + name)
    writer, name = generateModel(grid, model, *args, **kwargs)
    builder.load(name, writer.writer.getvalue(), "ellipticModel")
    writer.close()
    return importlib.import_module("dune.generated." + name)
