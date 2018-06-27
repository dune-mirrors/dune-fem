from __future__ import division, print_function, unicode_literals

from ufl import Form, Coefficient
from ufl.equation import Equation

from dune.common.compatibility import isString

from dune.source.cplusplus import Include, NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter
from dune.source.fem import fieldTensorType
from dune.ufl import GridFunction, DirichletBC
from dune.ufl.gatherderivatives import gatherDerivatives

from .ufl import compileUFL, fieldVectorType

def init(integrands, *args, **kwargs):
    coefficients = kwargs.pop('coefficients', dict())
    if len(args) == 1 and isinstance(args[0], dict):
        coefficients.update(args[0])
        args = []
    else:
        args = list(a for a in args if not isinstance(a,DirichletBC))
        dirichletBCs = list(a for a in args if isinstance(a,DirichletBC))

    coefficientNames = integrands._coefficientNames
    if len(args) > len(coefficientNames) + len(dirichletBCs):
        raise ValueError('Too many coefficients passed.')
    args += [None] * (len(coefficientNames) - len(args))

    for name, value in kwargs.items():
        try:
            i = coefficientNames[name]
        except KeyError:
            raise ValueError('No such coefficent: ' + name + '.')

        if args[i] is not None:
            raise ValueError('Coefficient already given as positional argument: ' + name + '.')
        args[i] = value

    for key, value in coefficients.items():
        if isinstance(key, Coefficient):
            try:
                i = integrands._renumbering[key]
            except AttributeError:
                raise ValueError('Cannot map UFL coefficients, because model was not generated from UFL form.')
            except KeyError:
                raise ValueError('No such coefficient: ' + str(key) + '.')
        elif isString(key):
            try:
                i = coefficientNames[key]
            except KeyError:
                raise ValueError('No such coefficent: ' + key + '.')
        else:
            raise ValueError('Expecting keys of coefficient map to be strings or  intances of ufl.Coefficient.')
        if args[i] is not None:
            raise ValueError('Coefficient already given as positional or keyword argument: ' + str(key) + '.')
        args[i] = value

    if hasattr(integrands, '_renumbering'):
        for c, i in integrands._renumbering.items():
            if isinstance(c, GridFunction):
                if args[i] is None:
                    args[i] = c.gf

    if any(arg is None for arg in args):
        missing = [name for name, i in coefficientNames.items() if args[i] is None]
        raise ValueError('Missing coefficients: ' + ', '.join(missing) + '.')

    integrands._init(*args)


def setConstant(integrands, index, value):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setConstant(index, value)


class Source(object):
    def __init__(self, gridType, gridIncludes, integrands,*args,
            tempVars=True, virtualize=True):
        self.gridType = gridType
        self.gridIncludes = gridIncludes
        self.integrands = integrands
        self.tempVars = tempVars
        self.virtualize = virtualize
        self.args = args

    def signature(self):
        return self.integrands.signature()

    def name(self):
        from dune.common.hashit import hashIt
        if self.virtualize:
            return 'integrands_' + self.signature() + '_' + hashIt(self.gridType)
        else:
            return 'integrands_nv' + self.signature() + '_' + hashIt(self.gridType)

    def valueTuples(self):
        if isinstance(self.integrands, Form):
            derivatives = gatherDerivatives(self.integrands)
            return ['std::tuple< ' + ', '.join(fieldTensorType(v.ufl_shape) for v in d) + ' >' for d in derivatives]
        else:
            return [self.integrands.rangeValueTuple, self.integrands.domainValueTuple]

    def __str__(self):
        if isinstance(self.integrands, Form):
            coefficients = set(self.integrands.coefficients())
            constants = [c for c in coefficients if c.is_cellwise_constant()]
            coefficients = sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())
            integrands = compileUFL(self.integrands, *self.args, constants=constants, coefficients=coefficients, tempVars=self.tempVars)
        else:
            integrands = self.integrands

        code = [Include('config.h')]
        code += [Include(i) for i in self.gridIncludes]

        code += integrands.includes()
        code.append(Include("dune/python/pybind11/pybind11.h"))
        code.append(Include("dune/python/pybind11/extensions.h"))
        code.append(Include("dune/fempy/py/grid/gridpart.hh"))

        if integrands._coefficients:
            if self.virtualize:
                code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
            else:
                for c in coefficients:
                    for i in c._includes:
                        code.append(Include(i))
        code.append(Include("dune/fempy/py/integrands.hh"))

        nameSpace = NameSpace('Integrands_' + integrands.signature())
        nameSpace.append(integrands.code())
        code.append(nameSpace)

        code.append(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + self.gridType + ' >'))
        if integrands._coefficients:
            if self.virtualize:
                coefficients = ['Dune::FemPy::VirtualizedGridFunction< GridPart, ' + c + ' >' for c in integrands._coefficients]
            else:
                coefficients = [c._typeName for c in coefficients]
        else:
            coefficients = []
        integrandsName = nameSpace.name + '::Integrands< ' + ', '.join(['GridPart'] + coefficients) + ' >'
        code.append(TypeAlias('Integrands', integrandsName))

        writer = SourceWriter()
        writer.emit(code);

        name = self.name()
        writer.openPythonModule(name)
        writer.emit('auto cls = Dune::Python::insertClass<Integrands>(module,"Integrands",Dune::Python::GenerateTypeName("'+integrandsName+'"), Dune::Python::IncludeFiles({"python/dune/generated/'+name+'.cc"})).first;')
        writer.emit('Dune::FemPy::registerIntegrands< Integrands >( module, cls );')
        if coefficients:
            coefficientNames = integrands.coefficientNames
            initArgs = ', '.join('const ' + t + ' &' + n for t, n in zip(coefficients, coefficientNames))
            keepAlive = ', '.join('pybind11::keep_alive< 1, ' + str(i+2) + ' >()' for i in range(len(coefficientNames)))
            writer.emit('cls.def( pybind11::init( [] ( ' + initArgs + ' ) { return new Integrands( ' + ', '.join(coefficientNames) +  ' ); } ), ' + keepAlive + ' );')
        else:
            writer.emit('cls.def( pybind11::init( [] () { return new Integrands(); } ) );')

        for t, n in zip(integrands.constantTypes, integrands.constantNames):
            te = "Integrands::" + t
            writer.emit('cls.def_property( "' + n + '", [] ( Integrands &self ) -> ' + te + ' { return self.' + n + '(); }, [] ( Integrands &self, const ' + te + ' &v ) { self.' + n + '() = v; } );')
        writer.emit('cls.def_property_readonly( "virtualized", [] ( Integrands& ) -> bool { return '+str(self.virtualize).lower()+';});')
        hasDirichletBC = 'true' if integrands.hasDirichletBoundary else 'false'
        writer.emit('cls.def_property_readonly( "hasDirichletBoundary", [] ( Integrands& ) -> bool { return '+hasDirichletBC+';});')

        writer.closePythonModule(name)

        source = writer.writer.getvalue()
        writer.close()
        return source


def load(grid, integrands, *args, renumbering=None, tempVars=True, virtualize=True):
    if isinstance(integrands, Equation):
        integrands = integrands.lhs - integrands.rhs

    source = Source(grid._typeName, grid._includes, integrands,*args,
            tempVars=tempVars,virtualize=virtualize)
    if isinstance(integrands, Form):
        coefficients = set(integrands.coefficients())
        numCoefficients = len(coefficients)
        if renumbering is None:
            renumbering = dict()
            renumbering.update((c, i) for i, c in enumerate(sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())))
            renumbering.update((c, i) for i, c in enumerate(c for c in coefficients if c.is_cellwise_constant()))
        coefficientNames = ['coefficient' + str(i) if n is None else n for i, n in enumerate(getattr(c, 'name', None) for c in coefficients if not c.is_cellwise_constant())]
    else:
        coefficientNames = integrands.coefficientNames

    from dune.generator import builder
    module = builder.load(source.name(), source, "integrands")
    rangeValueTuple, domainValueTuple = source.valueTuples()
    setattr(module.Integrands, "_domainValueType", domainValueTuple)
    setattr(module.Integrands, "_rangeValueType", rangeValueTuple)
    if not hasattr(module.Integrands, "_init"):
        setattr(module.Integrands, '_coefficientNames', {n: i for i, n in enumerate(coefficientNames)})
        module.Integrands._init = module.Integrands.__dict__['__init__']
        setattr(module.Integrands, '__init__', init)
        if renumbering is not None:
            setattr(module.Integrands, '_renumbering', renumbering)
            module.Integrands._setConstant = module.Integrands.__dict__['setConstant']
            setattr(module.Integrands, 'setConstant', setConstant)
    return module
