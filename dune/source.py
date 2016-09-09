from __future__ import print_function

import sys


# FileWriter
# ----------

class FileWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")

    def emit(self, src):
        print(src, file=self.file)

    def close(self):
        self.file.close()


# ListWriter
# ----------

class ListWriter:
    def __init__(self):
        self.lines = []

    def emit(self, src):
        self.lines.append(src)

    def close(self):
        pass



# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, writer):
        if self._isstring(writer):
            self.writer = FileWriter(writer)
        else:
            self.writer = writer
        self.blocks = []
        self.begin = True

    def close(self):
        if self.blocks:
            raise Exception("Open blocks left in source file.")
        self.writer.close()

    def emit(self, src):
        if src is None:
            return
        elif isinstance(src, (list, set, tuple)):
            for srcline in src:
                self.emit(srcline)
        elif self._isstring(src):
            src = src.rstrip()
            if src:
                self.writer.emit('  ' * len(self.blocks) + src)
            else:
                self.writer.emit('')
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

    def getName(self,block):
        for b in reversed(self.blocks):
            if b[0] == block:
                return b[1]
        return None

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

    def openClass(self, name, targs=None, bases=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('class ' + name)
        if bases:
            for i in range(0,len(bases)):
                prefix = '  : ' if i == 0 else '    '
                postfix = ',' if i+1 < len(bases) else ''
                self.emit(prefix + base + postfix)
        self.emit('{')
        self.pushBlock('class', name)

    def closeClass(self, name=None):
        self.popBlock('class', name)
        self.emit('};')

    def openStruct(self, name, targs=None, bases=None):
        self.emit(None if self.begin else ['','',''])
        self.emit(['// ' + name, '// ' + '-' * len(name), ''])
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        self.emit('struct ' + name)
        if bases:
            for i in range(0,len(bases)):
                prefix = '  : ' if i == 0 else '    '
                postfix = ',' if i+1 < len(bases) else ''
                self.emit(prefix + bases[i] + postfix)
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

    def openConstMethod(self, typedName, targs=None, args=None, implemented=True):
        self.emit(None if self.begin else '')
        if targs:
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' ) const')
        else:
            self.emit(typedName + ' () const')
        self.emit('{')
        self.pushBlock('const method', typedName)
        if not implemented:
            struct = self.getName("struct")
            if struct:
                fullMethodName = struct + "::" + typedName
            else:
                fullMethodName = typedName
            self.emit('DUNE_THROW(Dune::NotImplemented, "'+fullMethodName+'");')

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

    def openRangeBasedFor(self, typedName, container):
        self.emit('for( ' + typedName.strip() + ' : ' + container.strip() + ' )')
        self.emit('{')
        self.pushBlock('range-based for', typedName)

    def closeRangeBasedFor(self, typedName=None):
        self.popBlock('range-based for', typedName)
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

    @staticmethod
    def cpp_fields(field):
        if field=="complex":
            return "std::complex<double>"
        else:
            return "double"

class BaseModel:
    def __init__(self, dimRange, signature):
        self.dimRange = dimRange
        self.coefficients = []
        self.init = None
        self.vars = None
        self.signature = signature
        self.field = "double"

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
                    [('Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' >')\
                    for coefficient in self.coefficients if coefficient['constant']]) + ' >',\
                    'ConstantsTupleType;')
            sourceWriter.typedef('std::tuple< ' + ', '.join(\
                    [('Dune::Fem::FunctionSpace< double,' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    'dimDomain, ' +\
                    str(coefficient['dimRange']) + ' >')\
                    for coefficient in self.coefficients if not coefficient['constant']]) + ' >',\
                    'CoefficientFunctionSpaceTupleType')

            sourceWriter.typedef('typename std::tuple_element_t<i,ConstantsTupleType>', 'ConstantsRangeType', targs=['std::size_t i'])
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
        sourceWriter.typedef('typename std::tuple_element< i, ConstantsTupleType >::type', 'ConstantsType', targs=['std::size_t i'])

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

    def post(self, sourceWriter, name='Model', targs=[]):
        sourceWriter.openConstMethod('const ConstantsType< i > &constant', targs=['std::size_t i'])
        sourceWriter.emit('return std::get< i >( constants_ );')
        sourceWriter.closeConstMethod()
        sourceWriter.openMethod('ConstantsType< i > &constant', targs=['std::size_t i'])
        sourceWriter.emit('return std::get< i >( constants_ );')
        sourceWriter.closeMethod()

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
        sourceWriter.emit('mutable ConstantsTupleType constants_;')
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

        sourceWriter.openFunction('void setConstant', targs=['std::size_t i'], args=[modelClass + ' &model', 'pybind11::list l'])
        sourceWriter.emit('model.template constant< i >() = l.template cast< typename ' + modelClass + '::ConstantsType<i> >();')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('auto defSetConstant', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.typedef('std::function< void( ' + modelClass + ' &model, pybind11::list ) >', 'Dispatch')
        sourceWriter.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setConstant< i > )... }};')
        sourceWriter.emit('')
        sourceWriter.emit('return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::list l ) {')
        sourceWriter.emit('    std::size_t k = coeff.attr("number").template cast<int>();')
        # sourceWriter.emit('    if ( !coeff("is_piecewise_constant")(coeff).template cast<bool>() )')
        # sourceWriter.emit('      throw std::range_error( "Using setConstant for a Coefficient" );' )
        sourceWriter.emit('    if( k >= dispatch.size() )')
        sourceWriter.emit('      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );' )
        sourceWriter.emit('    dispatch[ k ]( model.impl(), l );')
        sourceWriter.emit('    return k;')
        sourceWriter.emit('  };')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('void setCoefficient', targs=['std::size_t i'], args=[ modelClass + ' &model', 'pybind11::handle o'])
        sourceWriter.emit('model.template coefficient< i >() = o.template cast< typename std::tuple_element< i, Coefficients >::type >().localFunction();')
        sourceWriter.closeFunction()

        sourceWriter.openFunction('auto defSetCoefficient', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
        sourceWriter.typedef('std::function< void( ' + modelClass + ' &model, pybind11::handle ) >', 'Dispatch')
        sourceWriter.emit('std::array< Dispatch, sizeof...( i ) > dispatch = {{ Dispatch( setCoefficient< i > )... }};')
        sourceWriter.emit('')
        sourceWriter.emit('return [ dispatch ] ( ' + wrapperClass + ' &model, pybind11::handle coeff, pybind11::handle o ) {')
        sourceWriter.emit('    std::size_t k = coeff.attr("number").template cast<int>();')
        # sourceWriter.emit('    if ( coeff.attr("is_piecewise_constant").template cast<bool>() )')
        # sourceWriter.emit('      throw std::range_error( "Using setCoefficient for a Constant" );' )
        sourceWriter.emit('    if( k >= dispatch.size() )')
        sourceWriter.emit('      throw std::range_error( "No such coefficient: "+std::to_string(k)+" >= "+std::to_string(dispatch.size()) );' )
        sourceWriter.emit('    dispatch[ k ]( model.impl(), o );')
        sourceWriter.emit('    return k;')
        sourceWriter.emit('  };')
        sourceWriter.closeFunction()
