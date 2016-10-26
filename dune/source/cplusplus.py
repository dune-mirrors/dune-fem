from __future__ import print_function, unicode_literals

import copy, io, sys


# FileWriter
# ----------

class FileWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")

    def emit(self, src):
        print(src, file=self.file)

    def close(self):
        self.file.close()

class StringWriter:
    def __init__(self):
        self.file = io.StringIO()

    def emit(self, src):
        print(src, file=self.file)

    def getvalue(self):
        return self.file.getvalue()

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




def stripEach(l):
    return None if l is None else [o.strip() for o in l]



# Block
# -----

class Block:
    def __init__(self):
        self.content = []

    def append(self, *objs):
        for obj in objs:
            if isinstance(obj, (list, set, tuple)):
                self.content += [o for o in obj]
            else:
                self.content.append(obj)



# NameSpace
# ---------

class NameSpace(Block):
    def __init__(self, name=None):
        Block.__init__(self)
        self.name = name



# Function
# --------

class Function(Block):
    def __init__(self, typedName, targs=None, args=None, code=None):
        Block.__init__(self)
        self.typedName = typedName
        self.targs = stripEach(targs)
        self.args = stripEach(args)
        if code is not None:
            self.append(code)



# Method
# ------

class Method(Block):
    def __init__(self, typedName, targs=None, args=None, code=None, static=False, const=False, volatile=False):
        Block.__init__(self)
        self.typedName = typedName
        self.static = static
        self.resetQualifiers(const=const, volatile=volatile)
        if static and (const or volatile):
            raise Exception('Cannot cv-qualify static method.')
        self.targs = stripEach(targs)
        self.args = stripEach(args)
        if code is not None:
            self.append(code)

    def variant(self, typedName, const=False, volatile=False):
        method = copy.copy(self)
        method.typedName = typedName
        method.resetQualifiers(const=const, volatile=volatile)
        return method

    def resetQualifiers(self, const=False, volatile=False):
        self.qualifiers=[]
        if const:
            self.qualifiers.append('const')
        if volatile:
            self.qualifiers.append('volatile')



# Constructor
# -----------

class Constructor(Block):
    def __init__(self, targs=None, args=None, init=None, code=None):
        Block.__init__(self)
        self.targs = stripEach(targs)
        self.args = stripEach(args)
        self.init = stripEach(init)
        if code is not None:
            self.append(code)



# TypeAlias
# ---------

class TypeAlias:
    def __init__(self, name, typeName, targs=None):
        self.name = name
        self.typeName = typeName
        self.targs = stripEach(targs)



# Variable
# --------

class Variable:
    def __init__(self, typedName, value=None, static=False, mutable=False):
        self.typedName = typedName
        self.value = value
        self.static = static
        self.mutable = mutable



# Class
# -----

class Class(Block):
    def __init__(self, name, targs=None, bases=None, final=False):
        Block.__init__(self)
        self.name = name
        self.targs = None
        self.targs = stripEach(targs)
        self.bases = stripEach(bases)
        self.final = final
        self.type = 'class'

    def addMethod(self, typedName, targs=None, args=None, static=False, const=False, volatile=False, implemented=True):
        method = Method(typedName, targs=targs, args=args, static=static, const=const, volatile=volatile)
        if not implemented:
            method.append('DUNE_THROW(Dune::NotImplemented, "' + name + '::' + typedName + '");')
        append(method)
        return method if implemented else None



# Struct
# ------

class Struct(Class):
    def __init__(self, name, targs=None, bases=None, final=False):
        Class.__init__(self, name, targs=targs, bases=bases, final=final)
        self.type = 'struct'



# AccessModifier
# --------------

class AccessModifier:
    def __init__(self, modifier):
        if modifier not in ['private', 'protected', 'public']:
            raise Exception(modifier + ' is not a valid access modifier.')
        self.modifier = modifier



# EnumClass
# ---------

class EnumClass:
    def __init__(self, name, values, base=None):
        self.name = name
        self.base = None if base is None else base.strip()
        self.values = stripEach(values)



# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, writer=None):
        if not writer:
            self.writer = StringWriter()
        elif self._isstring(writer):
            self.writer = FileWriter(writer)
        else:
            self.writer = writer
        self.blocks = []
        self.begin = True

    def close(self):
        if self.blocks:
            raise Exception("Open blocks left in source file.")
        self.writer.close()

    def emit(self, src, indent=0, context=None):
        if src is None:
            return
        elif isinstance(src, (list, set, tuple)):
            for s in src:
                self.emit(s, indent, context)
        elif isinstance(src, NameSpace):
            self.emit(None if self.begin else '')
            if src.name is not None:
                self.emit('namespace ' + src.name, indent)
            else:
                self.emit('namespace', indent)
            self.emit('{', indent)
            self.emit('', indent)
            self.emit(src.content, intent+1, src)
            self.emit('', indent)
            if src.name is not None:
                self.emit('} // namespace ' + src.name, indent)
            else:
                self.emit('} // anonymous namespace', indent)
        elif isinstance(src, Class):
            self.emit(None if self.begin else ['','',''], indent)
            self.emit(['// ' + src.name, '// ' + '-' * len(src.name), ''], indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            self.emit(src.type + ' ' + src.name + (' final' if src.final else ''))
            if src.bases:
                for i in range(0, len(src.bases)):
                    prefix = ': ' if i == 0 else '  '
                    postfix = ',' if i+1 < len(src.bases) else ''
                    self.emit(prefix + src.bases[i] + postfix, indent+1)
            self.emit('{', indent)
            self.emit(src.content, indent+1, src)
            self.emit('};', indent)
        elif isinstance(src, EnumClass):
            code = 'enum class ' + src.name
            if src.base is not None:
                code += ' : ' + src.base
            code += ' { ' + ', '.join(src.values) + ' };'
            self.emit(code, indent)
        elif isinstance(src, Function):
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = src.typedName + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            self.emit(signature, indent)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, Method):
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = ('static ' if src.static else '') + src.typedName + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            if src.qualifiers:
                signature += ' ' + ' '.join(src.qualifiers)
            self.emit(signature, indent)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, Constructor):
            if not isinstance(context, Class):
                raise Exception('Constructors can only occur in classes or structs')
            self.emit(None if self.begin else '', indent)
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
            signature = context.name + ' ('
            if src.args:
                signature += ' ' + ', '.join(src.args) + ' '
            signature += ')'
            self.emit(signature, indent)
            if src.init:
                for i in range(0,len(src.init)):
                    prefix = ': ' if i == 0 else '  '
                    postfix = ',' if i+1 < len(src.init) else ''
                    self.emit(prefix + src.init[i] + postfix, indent+1)
            if src.content:
              self.emit('{', indent)
              self.emit(src.content, indent+1, src)
              self.emit('}', indent)
            else:
              self.emit('{}', indent)
        elif isinstance(src, AccessModifier):
            if not isinstance(context, Class):
                raise Exception('Access modifiers can only occur in classes or structs')
            self.emit(None if self.begin else '')
            self.emit(src.modifier + ':', indent-1)
            self.begin = True
        elif isinstance(src, TypeAlias):
            if src.targs:
                self.emit('template< ' + ', '.join(src.targs) + ' >', indent)
                self.emit('using ' + src.name + ' = ' + src.typeName + ';', indent)
            else:
                self.emit('typedef ' + src.typeName + ' ' + src.name + ';', indent)
        elif isinstance(src, Variable):
            declaration = ('static ' if src.static else '') + ('mutable ' if src.mutable else '') + src.typedName
            if src.value is not None:
                declaration += ' = ' + src.value
            self.emit(declaration + ';', indent)
        elif self._isstring(src):
            src = src.rstrip()
            if src:
                self.writer.emit('  ' * (len(self.blocks) + indent) + src)
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
        self.emit(TypeAlias(typeAlias, typeName, targs=targs))

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
