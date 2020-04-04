from __future__ import division, print_function, unicode_literals

from ..common.utility import isString

from .expression import Variable

class BuiltInFunction:
    def __init__(self, header, cppType, name, namespace='std', targs=None, args=None):
        self.header = header
        self.cppType = cppType
        self.name = name
        self.namespace = namespace
        self.tarts = None if targs is None else [a.strip() for a in targs]
        self.args = [] if args is None else [a.strip() if isString(a) else a for a in args]

        if len(self.args) > 0:
            arg = self.args[len(self.args)-1]
            self.variadic = (isinstance(arg, Variable) and arg.cppType.endswith('...'))
        else:
            self.variadic = False

    def __call__(self, *args):
        if (len(args) != len(self.args)) and not self.variadic:
            raise Exception('Wrong number of Arguments: ' + str(len(args)) + ' (should be ' + str(len(self.args)) + ').')
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


and_ = BuiltInFunction('algorithm', 'X', 'std::logical_and', targs=['class X'], args=['const X &x', 'const X &y'])
max_ = BuiltInFunction('algorithm', 'X', 'max', targs=['class X'], args=['const X &x', 'const X &y'])
min_ = BuiltInFunction('algorithm', 'X', 'min', targs=['class X'], args=['const X &x', 'const X &y'])

abs_ = BuiltInFunction('cmath', 'X', 'abs', targs=['class X'], args=['const X &x'])
atan = BuiltInFunction('cmath', 'X', 'atan', targs=['class X'], args=['const X &x'])
atan2 = BuiltInFunction('cmath', 'X', 'atan2', targs=['class X'], args=['const X &x', 'const X &y'])
exp = BuiltInFunction('cmath', 'X', 'exp', targs=['class X'], args=['const X &x'])
cos = BuiltInFunction('cmath', 'X', 'cos', targs=['class X'], args=['const X &x'])
cosh = BuiltInFunction('cmath', 'X', 'cosh', targs=['class X'], args=['const X &x'])
log = BuiltInFunction('cmath', 'X', 'log', targs=['class X'], args=['const X &x'])
pow_ = BuiltInFunction('cmath', 'X', 'pow', targs=['class X'], args=['const X &x', 'const X &y'])
sin = BuiltInFunction('cmath', 'X', 'sin', targs=['class X'], args=['const X &x'])
sinh = BuiltInFunction('cmath', 'X', 'sinh', targs=['class X'], args=['const X &x'])
sqrt = BuiltInFunction('cmath', 'X', 'sqrt', targs=['class X'], args=['const X &x'])
tan = BuiltInFunction('cmath', 'X', 'tan', targs=['class X'], args=['const X &x'])
tanh = BuiltInFunction('cmath', 'X', 'tanh', targs=['class X'], args=['const X &x'])

def get(i):
    return BuiltInFunction('tuple', 'auto', 'get< ' + str(i) + ' >', targs=['class T'], args=['const T &arg'])

make_pair = BuiltInFunction('utility', 'std::pair< U, V >', 'make_pair', targs=['class U', 'class V'], args=['const U &left', 'const V &right'])

def make_index_sequence(n):
    return BuiltInFunction('utility', 'auto', 'make_index_sequence< ' + str(n) + ' >')

def make_shared(T):
    return BuiltInFunction('memory', 'std::shared_ptr< T >', 'make_shared< ' + T + ' >', targs=['class... Args'], args=[Variable('Args &&...', 'args')])

hybridForEach = BuiltInFunction('dune/common/hybridutilities.hh', 'void', 'forEach', namespace='Dune::Hybrid', targs=['class Range', 'class F'], args=['Range &&range', 'F &&f'])

coordinate = BuiltInFunction('dune/fem/common/coordinate.hh', 'X', 'coordinate', namespace='Dune::Fem', targs=['class X'], args=['const X &x'])

maxEdgeLength = BuiltInFunction('dune/fempy/geometry/edgelength.hh', 'typename Geometry::ctype', 'maxEdgeLength', namespace='Dune::FemPy', targs=['class Geometry'], args=['const Geometry &geometry'])
minEdgeLength = BuiltInFunction('dune/fempy/geometry/edgelength.hh', 'typename Geometry::ctype', 'minEdgeLength', namespace='Dune::FemPy', targs=['class Geometry'], args=['const Geometry &geometry'])
