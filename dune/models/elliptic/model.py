from __future__ import division, print_function, unicode_literals

from dune.source.builtin import get, hybridForEach, make_index_sequence, make_shared
from dune.source.cplusplus import UnformattedExpression
from dune.source.cplusplus import AccessModifier, Constructor, Declaration, Function, Method, NameSpace, Struct, TypeAlias, UnformattedBlock, Variable
from dune.source.cplusplus import assign, construct, dereference, lambda_, nullptr, return_
from dune.source.cplusplus import SourceWriter
from dune.source.fem import declareFunctionSpace


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
