from __future__ import division, print_function, unicode_literals

from dune.common.utility import isString

class Block:
    def __init__(self):
        self.content = []

    def append(self, *objs):
        for obj in objs:
            if isinstance(obj, (list, set, tuple)):
                self.content += [o for o in obj]
            elif obj is not None:
                self.content.append(obj)


class Statement:
    def __init__(self):
        pass


class UnformattedBlock():
    def __init__(self, *lines):
        self.lines = []
        self.append(*lines)

    def append(self, *lines):
        for line in lines:
            if isinstance(line, (list, tuple)):
                self.append(*list(line))
            elif isString(line):
                self.lines += line.split('\n')
            else:
                raise Exception('Only strings (or lists of them) can be appended to an UnformattedBlock')
