from __future__ import division, print_function, unicode_literals

import copy

from ..common.utility import isString

from .common import Block


class Function(Block):
    def __init__(self, cppType, name, targs=None, args=None, code=None):
        Block.__init__(self)
        self.cppType = cppType
        self.name = name
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() if isString(a) else a for a in args]
        if code is not None:
            self.append(code)


class Method(Block):
    def __init__(self, cppType, name, targs=None, args=None, code=None, static=False, const=False, inline=False, volatile=False):
        Block.__init__(self)
        self.cppType = cppType
        self.name = name
        self.static = static
        self.inline = inline
        self.resetQualifiers(const=const, volatile=volatile)
        if static and (const or volatile):
            raise Exception('Cannot cv-qualify static method.')
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() if isString(a) else a for a in args]
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


class Constructor(Block):
    def __init__(self, targs=None, args=None, init=None, code=None):
        Block.__init__(self)
        self.targs = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() if isString(a) else a for a in args]
        self.init = None if init is None else [a.strip() for a in init]
        if code is not None:
            self.append(code)
