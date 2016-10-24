from __future__ import print_function, unicode_literals

from .cplusplus import Method, SourceWriter, Variable

class BaseModel:
    def __init__(self, dimRange, signature):
        self.dimRange = dimRange
        self.coefficients = []
        self.init = None
        self.vars = None
        self.signature = signature
        self.field = "double"

    def getNumber(self, expr):
        e = [ ee for ee in self.coefficients if ee["name"] == str(expr) ]
        if len(e) > 1:
            raise KeyError('two coefficients provided with same name')
        return e[0]["number"]

    def pre(self, sourceWriter, name='Model', targs=[], bases=[]):
        sourceWriter.openStruct(name, targs=(['class GridPart'] + targs + ['class... Coefficients']),\
        bases=(bases))

        sourceWriter.typedef("GridPart", "GridPartType")
        sourceWriter.typedef(SourceWriter.cpp_fields(self.field), "RangeFieldType")

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

        if self.coefficients:
            sourceWriter.typedef('std::tuple< ' + ', '.join(\
                    [('std::shared_ptr<Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' >>')\
                    for coefficient in self.coefficients if coefficient['constant']]) + ' >',\
                    'ConstantsTupleType;')
            sourceWriter.typedef('std::tuple< ' + ', '.join(\
                    [('Dune::Fem::FunctionSpace< double,' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    'dimDomain, ' +\
                    str(coefficient['dimRange']) + ' >')\
                    for coefficient in self.coefficients if not coefficient['constant']]) + ' >',\
                    'CoefficientFunctionSpaceTupleType')

            sourceWriter.typedef('typename std::tuple_element_t<i,ConstantsTupleType>::element_type', 'ConstantsRangeType', targs=['std::size_t i'])
            sourceWriter.emit('static const std::size_t numCoefficients = std::tuple_size< CoefficientFunctionSpaceTupleType >::value;')
            sourceWriter.typedef('typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type', 'CoefficientFunctionSpaceType', targs=['std::size_t i'] )
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::RangeType', 'CoefficientRangeType', targs=['std::size_t i'])
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::JacobianRangeType', 'CoefficientJacobianRangeType', targs=['std::size_t i'])
            sourceWriter.typedef('typename CoefficientFunctionSpaceType< i >::HessianRangeType', 'CoefficientHessianRangeType', targs=['std::size_t i'])
        else:
            sourceWriter.emit('static const std::size_t numCoefficients = 0u;')
            sourceWriter.typedef('std::tuple<>', 'ConstantsTupleType')

        sourceWriter.emit('')
        sourceWriter.typedef('typename std::tuple_element< i, std::tuple< Coefficients... > >::type', 'CoefficientType', targs=['std::size_t i'])
        sourceWriter.typedef('typename std::tuple_element< i, ConstantsTupleType >::type::element_type', 'ConstantsType', targs=['std::size_t i'])

        sourceWriter.openMethod(name, args=[])
        sourceWriter.emit('constructConstants( std::make_index_sequence< std::tuple_size<ConstantsTupleType>::value >() );' )
        sourceWriter.closeMethod()

        init = Method('bool init', args=['const EntityType &entity'], const=True)
        init.append('entity_ = &entity;',
                    'initCoefficients( std::make_index_sequence< numCoefficients >() );',
                    self.init,
                    'return true;')

        entity = Method('const EntityType &entity', const=True)
        entity.append('return *entity_;')

        sourceWriter.emit([init, entity, Method('std::string name', const=True, code=['return "' + name + '";'])])

    def post(self, sourceWriter, name='Model', targs=[]):
        constant = Method('ConstantsType< i > &constant', targs=['std::size_t i'])
        constant.append('return *( std::get< i >( constants_ ) );')
        sourceWriter.emit([constant.variant('const ConstantsType< i > &constant', const=True), constant])

        coefficient = Method('CoefficientType< i > &coefficient', targs=['std::size_t i'])
        coefficient.append('return std::get< i >( coefficients_ );')
        sourceWriter.emit([coefficient.variant('const CoefficientType< i > &coefficient', const=True), coefficient])

        sourceWriter.section('private')

        initCoefficients = Method('void initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'], const=True)
        initCoefficients.append('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );')

        constructConstants = Method('void constructConstants', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        constructConstants.append('std::ignore = std::make_tuple( (std::get< i >( constants_ ) = std::make_shared<ConstantsType< i >>(), i)... );')

        sourceWriter.emit([initCoefficients, constructConstants])

        sourceWriter.emit('')
        sourceWriter.emit(Variable('const EntityType *entity_', 'nullptr', mutable=True))
        sourceWriter.emit(Variable('std::tuple< Coefficients... > coefficients_;', mutable=True))
        sourceWriter.emit(Variable('ConstantsTupleType constants_;', mutable=True))
        sourceWriter.emit(self.vars)
        sourceWriter.closeStruct(name)

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

        setCoefficient = Function('void setCoefficient', targs=['std::size_t i'], args=[ modelClass + ' &model', 'pybind11::handle o'])
        setCoefficient.append('model.template coefficient< i >() = o.template cast< typename std::tuple_element< i, Coefficients >::type >().localFunction();')
        sourceWriter.emit(setCoefficient)

        defSetCoefficient = Function('auto defSetCoefficient', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        defSetCoefficient.append(TypeAlias('Dispatch', 'std::function< void( ' + modelClass + ' &model, pybind11::handle ) >'),
                                 Variable('std::array< Dispatch, sizeof...( i ) > dispatch', '{{ Dispatch( setCoefficient< i > )... }}'),
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
