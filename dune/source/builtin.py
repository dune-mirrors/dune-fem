from __future__ import division, print_function, unicode_literals


class BuiltInFunction:
    def __init__(self, header, cppType, name, namespace='std', targs=None, args=None):
        self.header = header
        self.cppType = cppType
        self.name = name
        self.namespace = namespace
        self.tarts = None if targs is None else [a.strip() for a in targs]
        self.args = None if args is None else [a.strip() for a in args]

    def __call__(self, *args):
        if len(args) != len(self.args):
            raise Exception('Wrong number of Arguments: ' + len(args) + ' (should be ' + len(self.args) + ').')
        from .expression import Application, makeExpression
        return Application(self, args=[makeExpression(arg) for arg in args])

    def __eq__(self, other):
        if not isinstance(other, BuiltInFunction):
            return False
        if (other.namespace == self.namespace) and (other.name == self.name):
            assert (other.header == self.header)
            return True
        else:
            return False

    def __hash__(self):
        return hash(("BuildInFunction", self.namespace, self.cppType, self.name))

    def __str__(self):
        if self.namespace is None:
            return self.name
        else:
            return self.namespace + '::' + self.name

    def __repr__(self):
        return "built-in(" + str(self) + ")"


atan = BuiltInFunction('cmath', 'X', 'atan', targs=['class X'], args=['conat X &x'])
atan2 = BuiltInFunction('cmath', 'X', 'atan2', targs=['class X'], args=['const X &x', 'const X &y'])
cos = BuiltInFunction('cmath', 'X', 'cos', targs=['class X'], args=['const X &x'])
pow_ = BuiltInFunction('cmath', 'X', 'pow', targs=['class X'], args=['const X &x', 'const X &y'])
sin = BuiltInFunction('cmath', 'X', 'sin', targs=['class X'], args=['const X &x'])
tan = BuiltInFunction('cmath', 'X', 'tan', targs=['class X'], args=['const X &x'])

def get(i):
    return BuiltInFunction('tuple', 'auto', 'get< ' + str(i) + ' >', targs=['class T'], args=['const T &arg'])

make_pair = BuiltInFunction('utility', 'std::pair< U, V >', 'make_pair', targs=['class U', 'class V'], args=['const U &left', 'const V &right'])

coordinate = BuiltInFunction('dune/fem/common/coordinate.hh', 'X', 'coordinate', namespace='Dune::Fem', targs=['class X'], args=['const X &x'])
