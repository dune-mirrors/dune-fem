assert False, "not used?"

from __future__ import print_function, unicode_literals

from .cplusplus import AccessModifier, Constructor, Declaration, Function, Method, TypeAlias, UnformattedExpression, UnformattedBlock, Variable
from .cplusplus import SourceWriter
from .cplusplus import nullptr
from .fem import declareFunctionSpace

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

    def pre(self, name='Model'):
        code = [TypeAlias('GridPartType', 'GridPart'),
                TypeAlias("EntityType", "typename GridPart::template Codim< 0 >::EntityType"),
                TypeAlias("IntersectionType", "typename GridPart::IntersectionType")]
        code.append(declareFunctionSpace("typename GridPartType::ctype", SourceWriter.cpp_fields(self.field), UnformattedExpression("int", "GridPartType::dimensionworld"), self.dimRange))
        code.append(Declaration(Variable("const int", "dimLocal"), initializer=UnformattedExpression("int", "GridPartType::dimension"), static=True))

        if self.coefficients:
            code.append(TypeAlias('ConstantsTupleType', 'std::tuple< ' + ', '.join(\
                    [('std::shared_ptr<Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' >>')\
                    for coefficient in self.coefficients if coefficient['constant']]) + ' >'))
            code.append(TypeAlias('CoefficientFunctionSpaceTupleType', 'std::tuple< ' + ', '.join(\
                    [('Dune::Fem::FunctionSpace< double,' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    'dimDomain, ' +\
                    str(coefficient['dimRange']) + ' >')\
                    for coefficient in self.coefficients if not coefficient['constant']]) + ' >'))

            code.append(TypeAlias('ConstantsRangeType', 'typename std::tuple_element_t<i,ConstantsTupleType>::element_type', targs=['std::size_t i']))
            code.append(Declaration(Variable('const std::size_t', 'numCoefficients'), initializer=UnformattedExpression('std::size_t', 'std::tuple_size< CoefficientFunctionSpaceTupleType >::value'), static=True))
            code.append(TypeAlias('CoefficientFunctionSpaceType', 'typename std::tuple_element< i, CoefficientFunctionSpaceTupleType >::type', targs=['std::size_t i']))
            code.append(TypeAlias('CoefficientRangeType', 'typename CoefficientFunctionSpaceType< i >::RangeType', targs=['std::size_t i']))
            code.append(TypeAlias('CoefficientJacobianRangeType', 'typename CoefficientFunctionSpaceType< i >::JacobianRangeType', targs=['std::size_t i']))
            code.append(TypeAlias('CoefficientHessianRangeType', 'typename CoefficientFunctionSpaceType< i >::HessianRangeType', targs=['std::size_t i']))
        else:
            code.append(Declaration(Variable('const std::size_t', 'numCoefficients'), initializer=0, static=True))
            code.append(TypeAlias('ConstantsTupleType', 'std::tuple<>'))


        code.append(TypeAlias('CoefficientType', 'typename std::tuple_element< i, std::tuple< Coefficients... > >::type', targs=['std::size_t i']))
        code.append(TypeAlias('ConstantsType', 'typename std::tuple_element< i, ConstantsTupleType >::type::element_type', targs=['std::size_t i']))

        code.append(Constructor(code=UnformattedBlock('constructConstants( std::make_index_sequence< std::tuple_size<ConstantsTupleType>::value >() );')))

        init = Method('bool', 'init', args=['const EntityType &entity'], const=True)
        init.append('entity_ = &entity;',
                    'initCoefficients( std::make_index_sequence< numCoefficients >() );',
                    self.init,
                    'return true;')

        entity = Method('const EntityType &', 'entity', const=True)
        entity.append('return *entity_;')

        code.append([init, entity])
        code.append(Method('std::string', 'name', const=True, code=['return "' + name + '";']))
        return code

    def post(self):
        constant = Method('ConstantsType< i > &', 'constant', targs=['std::size_t i'])
        constant.append('return *( std::get< i >( constants_ ) );')

        code = [constant.variant('const ConstantsType< i > &constant', const=True), constant]

        coefficient = Method('CoefficientType< i > &', 'coefficient', targs=['std::size_t i'])
        coefficient.append('return std::get< i >( coefficients_ );')
        code +=[coefficient.variant('const CoefficientType< i > &coefficient', const=True), coefficient]

        code.append(AccessModifier('private'))

        initCoefficients = Method('void', 'initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'], const=True)
        initCoefficients.append('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).bind( entity() ), i)... );')

        constructConstants = Method('void', 'constructConstants', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        constructConstants.append('std::ignore = std::make_tuple( (std::get< i >( constants_ ) = std::make_shared<ConstantsType< i >>(), i)... );')

        code += [initCoefficients, constructConstants]

        code.append(Declaration(Variable('const EntityType *', 'entity_'), nullptr, mutable=True))
        code.append(Declaration(Variable('std::tuple< Coefficients... >', 'coefficients_;'), mutable=True))
        code.append(Declaration(Variable('ConstantsTupleType', 'constants_;'), mutable=True))
        code.append(self.vars)
        return code

    def export(self, sourceWriter, modelClass='Model', wrapperClass='ModelWrapper', constrArgs=(), constrKeepAlive=None):
        if self.coefficients:
            # sourceWriter.emit('cls.def( "setCoefficient", defSetCoefficient( std::make_index_sequence< std::tuple_size<Coefficients>::value >() ) );')
            sourceWriter.emit('cls.def( "setConstant", defSetConstant( std::make_index_sequence< std::tuple_size <typename '+ modelClass + '::ConstantsTupleType>::value >() ) );')
        sourceWriter.emit('')
        sourceWriter.emit('cls.def( pybind11::init( [] ('+\
                ' '.join( i[1]+' '+i[0]+',' for i in constrArgs) +\
                'const pybind11::dict &coeff) {')
        sourceWriter.emit('  auto instance = new ' + wrapperClass + '('+\
                ', '.join(i[0] for i in constrArgs) +\
                ');')
        if self.coefficients:
            sourceWriter.emit('  const int size = std::tuple_size<Coefficients>::value;')
            sourceWriter.emit('  auto dispatch = defSetCoefficient( std::make_index_sequence<size>() );' )
            sourceWriter.emit('  std::vector<bool> coeffSet(size,false);')
            sourceWriter.emit('  for (auto item : coeff) {')
            sourceWriter.emit('    int k = dispatch(*instance, item.first, item.second); ')
            sourceWriter.emit('    coeffSet[k] = true;')
            sourceWriter.emit('  }')
            sourceWriter.emit('  if ( !std::all_of(coeffSet.begin(),coeffSet.end(),[](bool v){return v;}) )')
            sourceWriter.emit('    throw pybind11::key_error("need to set all coefficients during construction");')
            sourceWriter.emit('  return instance;')
        sourceWriter.emit('  return instance;')
        sourceWriter.emit('  }),')
        if constrKeepAlive:
            sourceWriter.emit(constrKeepAlive + ',')
        sourceWriter.emit(''.join('pybind11::arg("' + i[0] + '"), ' for i in constrArgs) +\
                'pybind11::arg("coefficients")=pybind11::dict() );')

    def codeCoefficient(self, code, coefficients, constants):
        """find coefficients/constants in code string and do replacements
        """
        if coefficients:
            l = list(coefficients.items())
            l.sort()
            number = 0
            numCoeffs = [coef['number'] for coef in self.coefficients if coef['constant'] == False]
            if numCoeffs:
                number = max(numCoeffs) + 1
            for key,val in l:
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
