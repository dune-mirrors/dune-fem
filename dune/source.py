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
        if _self._isstring(writer):
            self.writer = FileWriter(fileName)
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
                self.emit(prefix + base + postfix)
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
            self.emit('template< ' + ', '.join([arg.strip() for arg in targs]) + ' >')
        if args:
            self.emit(typedName + ' ( ' + ', '.join([arg.strip() for arg in args]) + ' ) const')
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
