class Operator:
    def __init__(self, name, numArgs):
        self.name = name
        self.numArgs = numArgs


class PrefixUnaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 1)


class PostfixUnaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 1)


class BinaryOperator(Operator):
    def __init__(self, name):
        Operator.__init__(self, name, 2)


class BracketOperator(Operator):
    def __init__(self):
        Operator.__init__(self, '[]', 2)
