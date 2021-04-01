from __future__ import division, print_function, unicode_literals

import re

from dune.common.utility import isInteger

from ufl import replace
from ufl.log import UFLException
from ufl.core.expr import Expr

from dune.source.builtin import get, make_shared
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, dereference, lambda_, nullptr, return_, this
from dune.source.cplusplus import SourceWriter
from dune.source.fem import declareFunctionSpace
from dune.ufl.codegen import generateMethod

from ufl.differentiation import Grad

class EllipticModel:
    version = "v1_1"
    def __init__(self, dimDomain, dimRange, u, signature):
        assert isInteger(dimRange)
        self.dimDomain = dimDomain
        self.dimRange = dimRange
        self.trialFunction = u
        self.init = []
        self.vars = None
        self.signature = signature + EllipticModel.version
        self.field = "double"

        self._constants = []
        self._constantNames = {}
        self._parameterNames = {}
        self._coefficients = []

        self.arg_r = Variable("RRangeType &", "result")
        self.arg_dr = Variable("RJacobianRangeType &", "result")

        self.arg_x = Variable("const Point &", "x")
        self.arg_u = Variable("const DRangeType &", "u")
        self.arg_du = Variable("const DJacobianRangeType &", "du")
        self.arg_d2u = Variable("const DHessianRangeType &", "d2u")
        self.arg_ubar = Variable("const DRangeType &", "ubar")
        self.arg_dubar = Variable("const DJacobianRangeType &", "dubar")
        self.arg_d2ubar = Variable("const DHessianRangeType &", "d2ubar")

        self.arg_i = Variable("const IntersectionType &", "intersection")
        self.arg_bndId = Variable("int", "bndId")

        self.source = [assign(self.arg_r, construct("RRangeType", 0))]
        self.linSource = [assign(self.arg_r, construct("RRangeType", 0))]
        self.diffusiveFlux = [assign(self.arg_dr, construct("RJacobianRangeType", 0))]
        self.linDiffusiveFlux = [assign(self.arg_dr, construct("RJacobianRangeType", 0))]
        self.fluxDivergence = [assign(self.arg_r, construct("RRangeType", 0))]
        self.alpha = [assign(self.arg_r, construct("RRangeType", 0))]
        self.linAlpha = [assign(self.arg_r, construct("RRangeType", 0))]

        self.hasDirichletBoundary = False
        self.hasNeumanBoundary = False
        self.isDirichletIntersection = [return_(False)]
        self.dirichlet = [assign(self.arg_r, construct("RRangeType", 0))]
        self.symmetric = False

        self.baseName = "elliptic"
        self.modelWrapper = "DiffusionModelWrapper< Model >"

    def predefineCoefficients(self,predefined,x):
        for coefficient, idx in self.coefficients.items():
            for derivative in self.coefficient(idx, x):
                predefined[coefficient] = derivative
                coefficient = Grad(coefficient)
        predefined.update({c: self.constant(i) for c, i in self.constants.items()})

    def addCoefficient(self, dimRange, typeName, name=None, field="double"):
        idx = len(self._coefficients)
        self._coefficients.append({'typeName':typeName, 'dimRange': dimRange, 'name': name, 'field': field})
        return idx

    def addConstant(self, cppType, name=None, parameter=None):
        idx = len(self._constants)
        self._constants.append(cppType)
        if name is not None:
            self._constantNames[name] = idx
        if parameter is not None:
            self._parameterNames[parameter] = idx
        return idx

    def cppIdentifier(self,name,base,idx):
        if re.match('^[a-zA-Z_][a-zA-Z0-9_]*$', name) is None:
            return base+str(idx)
        else:
            return name
    def cppTypeIdentifier(self,name,base,idx):
        ret = self.cppIdentifier(name,base,idx)
        return ret[0].upper() + ret[1:]
    def cppVarIdentifier(self,name,base,idx):
        ret = self.cppIdentifier(name,base,idx)
        return ret[0].lower() + ret[1:]

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def constant(self, idx):
        return UnformattedExpression(self._constants[idx], 'constant< ' + str(idx) + ' >()')

    def coefficient(self, idx, x):
        coefficient = []
        for t, n in (('RangeType', 'evaluate'), ('JacobianRangeType', 'jacobian'), ('HessianRangeType', 'hessian')):
            result = Variable('typename std::tuple_element_t< ' + str(idx) + ', CoefficientFunctionSpaceTupleType >::' + t, 'result')
            code = [Declaration(result),
                    UnformattedExpression('void', 'std::get< ' + str(idx) + ' >( coefficients_ ).' + n + '( x, ' + result.name + ' )', uses=[result]),
                    return_(result)]
            coefficient += [lambda_(capture=[this], args=['auto x'], code=code)(x)]
        return coefficient

    @property
    def hasCoefficients(self):
        return bool(self._coefficients)

    @property
    def hasConstants(self):
        return bool(self._constants)

    def code(self, name=None, targs=None):
        if targs is None:
            targs = []
        if name is None:
            name = 'Model'
        constants_ = Variable('std::tuple< ' + ', '.join('std::shared_ptr< ' + c  + ' >' for c in self._constants) + ' >', 'constants_')
        # coefficients_ = Variable('std::tuple< ' + ', '.join(c['name'] if c['name'] is not None else 'Coefficient' + str(i) for i, c in enumerate(self._coefficients)) + ' >', 'coefficients_')
        coefficients_ = Variable('std::tuple< ' + ', '.join(\
                'Dune::Fem::ConstLocalFunction<' + self.cppTypeIdentifier(c['name'],"coefficient",i) + '> ' for i, c in enumerate(self._coefficients)) + ' >', 'coefficients_')
        entity_ = Variable('const EntityType *', 'entity_')

        # code = Struct(name, targs=(['class GridPart'] + ['class ' + c['name'] if c['name'] is not None else 'class Coefficient' + str(i) for i, c in enumerate(self._coefficients)] + targs))
        code = Struct(name, targs=(['class GridPart'] + ['class ' + self.cppTypeIdentifier(c['name'],"coefficient",i) for i, c in enumerate(self._coefficients)] + targs))

        code.append(TypeAlias("GridPartType", "GridPart"))
        code.append(TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"))
        code.append(TypeAlias("IntersectionType", "typename GridPart::IntersectionType"))

        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), UnformattedExpression("int", "GridPartType::dimensionworld"), self.dimDomain,
            name="DFunctionSpaceType",prefix="D",
            dimDomainName="dimDomain", dimRangeName="dimD"
            ))
        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), UnformattedExpression("int", "GridPartType::dimensionworld"), self.dimRange,
            name="RFunctionSpaceType",prefix="R",
            dimDomainName=None, dimRangeName="dimR"
            ))
        code.append(Declaration(Variable("const int", "dimLocal"), initializer=UnformattedExpression("int", "GridPartType::dimension"), static=True))

        if self.hasConstants:
            code.append(TypeAlias("ConstantType", "typename std::tuple_element_t< i, " + constants_.cppType + " >::element_type", targs=["std::size_t i"]))
            code.append(Declaration(Variable("const std::size_t", "numConstants"), initializer=len(self._constants), static=True))

        if self.hasCoefficients:
            code.append(TypeAlias('CoefficientType', 'std::tuple_element_t< i, ' + coefficients_.cppType + ' >', targs=['std::size_t i']))
            # coefficientSpaces = ["Dune::Fem::FunctionSpace< DomainFieldType, " + SourceWriter.cpp_fields(c['field']) + ", dimDomain, " + str(c['dimRange']) + " >" for c in self._coefficients]
            coefficientSpaces = ["typename CoefficientType<"+str(i)+">::FunctionSpaceType" for i,c in enumerate(self._coefficients)]
            code.append(TypeAlias("CoefficientFunctionSpaceTupleType", "std::tuple< " + ", ".join(coefficientSpaces) + " >"))

        arg_param = Variable("const Dune::Fem::ParameterReader &", "parameter")
        args = [Declaration(arg_param, initializer=UnformattedExpression('const ParameterReader &', 'Dune::Fem::Parameter::container()'))]
        init = None
        if self.hasCoefficients:
            # args = [Variable("const " + c['name'] if c['name'] is not None else "const Coefficient" + str(i) + " &", "coefficient" + str(i)) for i, c in enumerate(self._coefficients)] + args
            args = [Variable("const " + self.cppTypeIdentifier(c['name'],"coefficient",i) + " &", "coefficient" + str(i)) for i, c in enumerate(self._coefficients)] + args
            init = ["coefficients_(" + ",".\
                  join("CoefficientType<"+str(i)+">"\
                        +"(coefficient" + str(i)+")" for i, c in enumerate(self._coefficients)) + " )"]
        constructor = Constructor(args=args, init=init)
        constructor.append([assign(get(str(i))(constants_), make_shared(c)()) for i, c in enumerate(self._constants)])
        for name, idx in self._parameterNames.items():
            constructor.append(assign(dereference(get(idx)(constants_)), UnformattedExpression("auto", arg_param.name + '.getValue< ' + self._constants[idx] + ' >( "' + name + '" )', uses=[arg_param])))
        code.append(constructor)

        init = ['entity_ = &entity;']
        init += ['std::get< ' + str(i) + ' >( ' + coefficients_.name + ').bind( entity );' for i, c in enumerate(self._coefficients)]
        init = [UnformattedBlock(init)] + self.init + [return_(True)]
        code.append(Method('bool', 'init', args=['const EntityType &entity'], code=init, const=True))

        code.append(Method('const EntityType &', 'entity', code=return_(dereference(entity_)), const=True))
        code.append(Method('std::string', 'name', const=True, code=return_(UnformattedExpression('const char *', '"' + name + '"'))))

        code.append(TypeAlias("BoundaryIdProviderType", "Dune::Fem::BoundaryIdProvider< typename GridPartType::GridType >"))
        code.append(Declaration(Variable("const bool", "symmetric"), initializer=self.symmetric, static=True))

        code.append(Method('void', 'source', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.source, const=True))
        code.append(Method('void', 'linSource', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_r], code=self.linSource, const=True))

        code.append(Method('void', 'diffusiveFlux', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.diffusiveFlux, const=True))
        code.append(Method('void', 'linDiffusiveFlux', targs=['class Point'], args=[self.arg_ubar, self.arg_dubar, self.arg_x, self.arg_u, self.arg_du, self.arg_dr], code=self.linDiffusiveFlux, const=True))

        code.append(Method('void', 'fluxDivergence', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_du, self.arg_d2u, self.arg_r], code=self.fluxDivergence, const=True))

        code.append(Method('void', 'alpha', targs=['class Point'], args=[self.arg_x, self.arg_u, self.arg_r], code=self.alpha, const=True))
        code.append(Method('void', 'linAlpha', targs=['class Point'], args=[self.arg_ubar, self.arg_x, self.arg_u, self.arg_r], code=self.linAlpha, const=True))

        code.append(Method('bool', 'hasNeumanBoundary', const=True, code=return_(self.hasNeumanBoundary)))

        code.append(TypeAlias("DirichletComponentType","std::array<int,"+str(self.dimRange)+">"))
        code.append(Method('bool', 'hasDirichletBoundary', const=True, code=return_(self.hasDirichletBoundary)))
        code.append(Method('bool', 'isDirichletIntersection', args=[self.arg_i, 'DirichletComponentType &dirichletComponent'], code=self.isDirichletIntersection, const=True))
        code.append(Method('void', 'dirichlet', targs=['class Point'], args=[self.arg_bndId, self.arg_x, self.arg_r], code=self.dirichlet, const=True))

        if self.hasConstants:
            code.append(Method("const ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_))), const=True))
            code.append(Method("ConstantType< i > &", "constant", targs=["std::size_t i"], code=return_(dereference(get("i")(constants_)))))

        if self.hasCoefficients:
            code.append(Method("const CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_)), const=True))
            code.append(Method("CoefficientType< i > &", "coefficient", targs=["std::size_t i"], code=return_(get("i")(coefficients_))))

        for n, i in self._constantNames.items():
            t = self._constants[i]
            code.append(Method('const ' + t + ' &', n, code=return_(dereference(get(i)(constants_))), const=True))
            code.append(Method(t + ' &', n, code=return_(dereference(get(i)(constants_)))))

        code.append(AccessModifier("private"))
        code.append(Declaration(entity_, nullptr, mutable=True))
        if self.hasConstants:
            code.append(Declaration(constants_, mutable=True))
        if self.hasCoefficients:
            code.append(Declaration(coefficients_, mutable=True))
        return code

    #def write(self, sourceWriter, name='Model', targs=[]):
    #    sourceWriter.emit(self.code(name=name, targs=targs))

    def exportSetConstant(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper'):
        sourceWriter.emit(TypeAlias('ModelType', modelClass))
        sourceWriter.openFunction('std::size_t renumberConstants', args=['pybind11::handle &obj'])
        sourceWriter.emit('std::string id = pybind11::str( obj );')
        sourceWriter.emit('if( pybind11::hasattr(obj,"name") ) id = pybind11::str(obj.attr("name"));')
        for name, number in self._constantNames.items():
            sourceWriter.emit('if (id == "' + name + '") return ' + str(number) + ';')
        sourceWriter.emit('throw pybind11::value_error("coefficient \'" + id + "\' has not been registered");')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('void setConstant', targs=['std::size_t i'], args=['ModelType &model', 'pybind11::handle value'])
        sourceWriter.emit('model.template constant< i >() = value.template cast< typename ModelType::ConstantType< i > >();')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('auto DUNE_PRIVATE defSetConstant', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.emit(TypeAlias('Dispatch', 'std::function< void( ModelType &model, pybind11::handle ) >'))
        sourceWriter.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setConstant< i > )... }};')
        sourceWriter.emit('')
        sourceWriter.emit('return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::handle value ) {')
        sourceWriter.emit('    std::size_t k = renumberConstants( coeff );')
        sourceWriter.emit('    if( k >= dispatch.size() )')
        sourceWriter.emit('      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );' )
        sourceWriter.emit('    dispatch[ k ]( model, value );')
        sourceWriter.emit('    return k;')
        sourceWriter.emit('  };')
        sourceWriter.closeFunction()

    def export(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper',nameSpace=''):
        if self.hasConstants:
            sourceWriter.emit('cls.def( "setConstant",'+nameSpace+'::defSetConstant( std::make_index_sequence< ' + modelClass + '::numConstants >() ) );')
        coefficients = [('Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< ' + SourceWriter.cpp_fields(c['field']) + ', ' + str(c['dimRange']) + ' > >')
                        if not c['typeName'].startswith("Dune::Python::SimpleGridFunction") \
                        else c['typeName'] \
                for c in self._coefficients]
        sourceWriter.emit('')
        # TODO
        sourceWriter.emit('cls.def( pybind11::init( [] ( ' + ', '.join( [] + ['const ' + c + ' &coefficient' + str(i) for i, c in enumerate(coefficients)]) + ' ) {')
        if self.hasCoefficients:
            sourceWriter.emit('  return new  ' + wrapperClass + '( ' + ', '.join('coefficient' + str(i) for i, c in enumerate(coefficients)) + ' );')
        else:
            sourceWriter.emit('  return new  ' + wrapperClass + '();')
        #if self.coefficients:
        #    sourceWriter.emit('  const int size = std::tuple_size<Coefficients>::value;')
        #    sourceWriter.emit('  auto dispatch = defSetCoefficient( std::make_index_sequence<size>() );' )
        #    sourceWriter.emit('  std::vector<bool> coeffSet(size,false);')
        #    sourceWriter.emit('  for (auto item : coeff) {')
        #    sourceWriter.emit('    int k = dispatch(instance, item.first, item.second); ')
        #    sourceWriter.emit('    coeffSet[k] = true;')
        #    sourceWriter.emit('  }')
        #    sourceWriter.emit('  if ( !std::all_of(coeffSet.begin(),coeffSet.end(),[](bool v){return v;}) )')
        #    sourceWriter.emit('    throw pybind11::key_error("need to set all coefficients during construction");')
        if self.hasCoefficients:
            sourceWriter.emit('  }), ' + ', '.join('pybind11::keep_alive< 1, ' + str(i) + ' >()' for i, c in enumerate(coefficients, start=2)) + ' );')
        else:
            sourceWriter.emit('  } ) );')

    def codeCoefficient(self, code, coefficients, constants):
        """find coefficients/constants in code string and do replacements
        """
        for name, value in coefficients.items():
            if not any(name == c['name'] for c in self._coefficients):
                self.addCoefficient(value.dimRange,value._typeName, name)
        for name, dimRange in constants.items():
            if name not in self._constantNames:
                self.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >', name)

        if '@const:' in code:
            codeCst = code.split('@const:')
            import itertools
            for name in set([''.join(itertools.takewhile(str.isalpha, str(c))) for c in codeCst[1:]]):
                if name not in self._constantNames:
                    cname = '@const:' + name
                    afterName = code.split(cname)[1:]
                    if afterName[0][0] == '[':
                        beforeText = [an.split(']')[0].split('[')[1] for an in afterName]
                        dimRange = max( [int(bt) for bt in beforeText] ) + 1
                    else:
                        dimRange = 1
                    self.addConstant('Dune::FieldVector< double, ' + str(dimRange) + ' >', name)

        for i, c in enumerate(self._coefficients):
            jacname = '@jac:' + c['name']
            if jacname in code:
                varname = 'dc' + str(i)
                code = code.replace(jacname, varname)
                decl = 'CoefficientJacobianRangeType< ' + str(i) + ' > ' + varname + ';'
                if not decl in code:
                    code = decl + '\ncoefficient< ' + str(i) + ' >().jacobian( x, ' + varname + ' );' + code
            gfname = '@gf:' + c['name']
            if gfname in code:
                varname = 'c' + str(i)
                code = code.replace(gfname, varname)
                decl = 'CoefficientRangeType< ' + str(i) + ' > c' + str(i) + ';'
                if not decl in code:
                    code = decl + '\ncoefficient< ' + str(i) + ' >().evaluate( x, ' + varname + ' );' + code

        for name, i in self._constantNames.items():
            cname = '@const:' + name
            if cname in code:
                varname = 'cc' + str(i)
                code = code.replace(cname, varname)
                init = 'const ' + self._constants[i] + ' &' + varname + ' = constant< ' + str(i) + ' >();'
                if not init in code:
                    code = init + code

    def appendCode(self, key, code, **kwargs):
        function = getattr(self, key)
        coef = kwargs.pop("coefficients", {})
        const = kwargs.pop("constants", {})
        function.append(UnformattedBlock(self.codeCoefficient(code, coef, const)))
        setattr(self, key, function)

    def generateMethod(self, code, expr, *args, **kwargs):
        if isinstance(expr, Expr):
            try:
                expr = replace(expr, self._replaceCoeff)
            except UFLException:
                pass
        return generateMethod(code,expr,*args,**kwargs)
