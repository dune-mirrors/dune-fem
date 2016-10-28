from __future__ import division, print_function, unicode_literals


class Block:
    def __init__(self):
        self.content = []

    def append(self, *objs):
        for obj in objs:
            if isinstance(obj, (list, set, tuple)):
                self.content += [o for o in obj]
            else:
                self.content.append(obj)


class Statement:
    def __init__(self):
        pass
