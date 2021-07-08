from __future__ import division, print_function, unicode_literals

import io, sys

from ..common.utility import isString

from .builtin import *
from .common import *
from .expression import *
from .function import *
from .operator import *

from .formatter.expression import formatExpression


# FileWriter
# ----------

class FileWriter:
    def __init__(self, fileName):
        self.file = open(fileName, "wt")

    def emit(self, src):
        print(src, file=self.file)

    def close(self):
        self.file.close()



# StringWriter
# ------------

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



class Include:
    def __init__(self, fileName):
        self.fileName = fileName
    def __lt__(self,other):
        return self.fileName < other.fileName


# NameSpace
# ---------

class NameSpace(Block):
    def __init__(self, name=None, code=None):
        Block.__init__(self)
        self.name = name
        if code is not None:
            self.append(code)



# TypeAlias
# ---------

class TypeAlias:
    def __init__(self, name, typeName, targs=None):
        self.name = name
        self.typeName = typeName
        self.targs = None if targs is None else [a.strip() for a in targs]



# Declaration
# -----------

class Declaration:
    def __init__(self, obj, initializer=None, static=False, mutable=False, constexpr=False):
        if not isinstance(obj, Variable):
            raise Exception('Only variables can be declared for now.')
        self.obj = obj
        self.initializer = None if initializer is None else makeExpression(initializer)
        self.static = static
        self.constexpr = constexpr
        self.mutable = mutable

    def __repr__(self):
        s = "declare " + self.obj.name + " : " + self.obj.cppType
        if self.initializer is not None:
            s += " = " + repr(self.initializer)
        return s



# Using
# -----

class Using:
    def __init__(self, obj):
        # if not isinstance(obj, BuiltInFunction):
        #     raise Exception('Only built-in functions can be used for now')
        self.obj = obj

    def __eq__(self, other):
        return self.obj == other.obj

    def __hash__(self):
        return hash((Using, self.obj))

    def __str__(self):
        return "using " + str(self.obj) + ";"

    def __repr__(self):
        return "using " + repr(self.obj)

    def __lt__(self, other):
        return str(self.obj) < str(other.obj)


# Class
# -----

class Class(Block):
    def __init__(self, name, targs=None, bases=None, final=False):
        Block.__init__(self)
        self.name = name
        self.targs = None
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.bases = None if bases is None else [a.strip() for a in bases]
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
        self.values = [v.strip() for v in values]



# SwitchStatement
# ---------------

class SwitchStatement(Statement):
    def __init__(self, var, branches=None, default=None):
        self.var = var
        self.branches = dict()
        if branches is not None:
            for case, code in branch:
                self.append(case, code)
        self.default = Block()
        self.default.append(default)

    def append(self, case, code):
        if isinstance(case,bool):
            case = "true" if case else "false"
        block = self.branches.setdefault(case, Block())
        block.append(code)

# IfStatement
# ---------------

class IfStatement(Statement):
    def __init__(self, condition, thenBranch, elseBranch=None, constexpr=False):
        self.condition = condition
        self.thenBranch = Block()
        self.thenBranch.append(thenBranch)
        if elseBranch is not None:
            self.elseBranch = Block()
            self.elseBranch.append(elseBranch)
        else:
            self.elseBranch = None
        self.constexpr = constexpr

# ReturnStatement
# ---------------

class ReturnStatement(Statement):
    def __init__(self, expr=None):
        Statement.__init__(self)
        self.expression = None if expr is None else makeExpression(expr)



# short hand notation
# -------------------

def return_(expr=None):
    return ReturnStatement(expr)



# SourceWriter
# ------------

class SourceWriter:
    def __init__(self, writer=None):
        if not writer:
            self.writer = StringWriter()
        elif isString(writer):
            self.writer = FileWriter(writer)
        else:
            self.writer = writer
        self.blocks = []
        self.begin = True
        self.context = None

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
        elif isinstance(src, Include):
            self.emit('#include <' + src.fileName + '>')
        elif isinstance(src, NameSpace):
            self.emit(None if self.begin else '')
            if src.name is not None:
                self.emit('namespace ' + src.name, indent)
            else:
                self.emit('namespace', indent)
            self.emit('{', indent)
            self.emit('', indent)
            self.emit(src.content, indent+1, src)
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
            self.emit(src.type + ' ' + src.name + (' final' if src.final else ''), indent)
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
            signature = self.typedName(src) + ' ('
            if src.args:
                signature += ' ' + ', '.join(self.formatArgument(arg) for arg in src.args) + ' '
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
            signature = ('static ' if src.static else '') + ('inline ' if src.inline else '') + self.typedName(src) + ' ('
            if src.args:
                signature += ' ' + ', '.join(self.formatArgument(arg) for arg in src.args) + ' '
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
                signature += ' ' + ', '.join(self.formatArgument(arg) for arg in src.args) + ' '
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
        elif isinstance(src, Declaration):
            if not isinstance(src.obj, Variable):
                raise Exception('Only variables can be declared for now.')
            declaration = ('static ' if src.static else '') +\
                          ('constexpr ' if src.constexpr else '') +\
                          ('mutable ' if src.mutable else '') + self.typedName(src.obj)
            if src.initializer is not None:
                expr = formatExpression(src.initializer)
                expr[len(expr)-1] += ';'
                self.emit(declaration + ' = ' + expr[0], indent)
                for e in expr[1:]:
                    if isinstance(e, tuple):
                        self.emit(e, indent+2, Function('auto', '<lambda>'))
                    else:
                        self.emit(e, indent+1)
            else:
                self.emit(declaration + ';', indent)
        elif isinstance(src, Using):
            self.emit(str(src), indent)
        elif isinstance(src, Statement):
            if not isinstance(context, (Constructor, Function, Method)) and\
               not self.context in ["PythonModule"]:
                raise Exception('Statements can only occur in constructors, functions and methods')
            if isinstance(src, Expression):
                expr = formatExpression(src)
                expr[len(expr)-1] += ';'
                self.emit(expr[0], indent)
                for e in expr[1:]:
                    if isinstance(e, tuple):
                        self.emit(e, indent+2, Function('auto', '<lambda>'))
                    else:
                        self.emit(e, indent+1)
            elif isinstance(src, ReturnStatement):
                if src.expression is not None:
                    expr = formatExpression(src.expression)
                    expr[len(expr)-1] += ';'
                    self.emit('return ' + expr[0], indent)
                    for e in expr[1:]:
                        if isinstance(e, tuple):
                            self.emit(e, indent+2, Function('auto', '<lambda>'))
                        else:
                            self.emit(e, indent+1)
                else:
                    self.emit('return;', indent)
            elif isinstance(src, SwitchStatement):
                self.emit('switch( ' + src.var.name + ' )', indent)
                self.emit('{', indent)
                for case, code in src.branches.items():
                    self.emit('case ' + str(case) + ':', indent)
                    if code.content:
                        self.emit('{', indent+1)
                        self.emit(code.content, indent+2, context)
                        self.emit('}', indent+1)
                    self.emit('break;', indent+1)
                if src.default.content:
                    self.emit('default:', indent)
                    self.emit('{', indent+1)
                    self.emit(src.default.content, indent+2, context)
                    self.emit('}', indent+1)
                self.emit('}', indent)
            elif isinstance(src, IfStatement):
                self.emit(('if constexpr' if src.constexpr else 'if') +
                          '( ' + src.condition + ' )', indent)
                self.emit(src.thenBranch,indent)
                if src.elseBranch is not None:
                    self.emit('else', indent)
                    self.emit(src.elseBranch,indent)
            else:
                raise Exception('Unknown statement type.')
        elif isinstance(src, Block):
             self.emit('{', indent)
             self.emit(src.content, indent+2, context)
             self.emit('}', indent)
        elif isinstance(src, UnformattedBlock):
            if not isinstance(context, (Constructor, Function, Method)):
                raise Exception('UnformattedBlocks can only occur in constructors, functions and methods')
            if src.lines:
                self.emit('{', indent)
                for line in src.lines:
                    self.emit(line, indent+1, src)
                self.emit('}', indent)
            else:
                self.emit('{}', indent)
        elif isString(src):
            src = src.rstrip()
            if src:
                self.writer.emit('  ' * (len(self.blocks) + indent) + src)
            else:
                self.writer.emit('')
            self.begin = False
        else:
            raise Exception("Unable to print " + repr(src) + ".")

    def typedName(self, obj):
        if obj.cppType is None:
            raise Exception('object ' + obj.name + '  does not have a type.')
        if obj.cppType.endswith('&') or obj.cppType.endswith('*'):
            return obj.cppType + obj.name
        else:
            return obj.cppType + ' ' + obj.name

    def formatArgument(self, arg):
        if isinstance(arg, Variable):
            return self.typedName(arg)
        elif isinstance(arg, Declaration):
            if arg.static:
                raise Exception("Arguments cannot be static.")
            if arg.mutable:
                raise Exception("Argumenta cannot be mutable.")
            if arg.constexpr:
                raise Exception("Argumenta cannot be constexpr.")
            declaration = self.typedName(arg.obj)
            if arg.initializer is not None:
                declaration += ' = ' + ' '.join(formatExpression(arg.initializer))
            return declaration
        elif isString(arg):
            return arg
        else:
            raise Exception("Unable to print argument " + repr(arg) + ".")

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

    def openPythonModule(self, moduleName):
        self.emit(None if self.begin else '')
        self.emit('PYBIND11_MODULE( ' + moduleName.strip() + ', module )')
        self.emit('{')
        self.pushBlock('pybind11 module', moduleName)
        self.context = "PythonModule"

    def closePythonModule(self, moduleName=None):
        self.popBlock('pybind11 module', moduleName)
        self.emit('}')
        self.context = None

    @staticmethod
    def cpp_fields(field):
        if field=="complex":
            return "std::complex<double>"
        else:
            return "double"
