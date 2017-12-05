from __future__ import division, print_function, unicode_literals

from ufl import Form
from ufl.equation import Equation

from dune.source.cplusplus import Include, NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter
from dune.source.fem import fieldTensorType
from dune.ufl.gatherderivatives import gatherDerivatives

from .ufl import compileUFL, fieldVectorType

def setConstant(integrands, index, value):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setConstant(index, value)


def setCoefficient(integrands, index, coefficient):
    try:
        index = integrands._renumbering[index]
    except KeyError:
        pass
    integrands._setCoefficient(index, coefficient)


class Source(object):
    def __init__(self, gridType, gridIncludes, integrands, tempVars=True):
        self.gridType = gridType
        self.gridIncludes = gridIncludes
        self.integrands = integrands
        self.tempVars = tempVars

    def signature(self):
        return self.integrands.signature()

    def name(self):
        from dune.common.hashit import hashIt
        return 'integrands_' + self.signature() + '_' + hashIt(self.gridType)


    def valueTuples(self):
        if isinstance(self.integrands, Form):
            derivatives = gatherDerivatives(self.integrands)
            return ['std::tuple< ' + ', '.join(fieldTensorType(v.ufl_shape) for v in d) + ' >' for d in derivatives]
        else:
            return [integrands.rangeValueType, integrands.domainValueTuple]

    def __str__(self):
        if isinstance(self.integrands, Form):
            coefficients = set(self.integrands.coefficients())
            constants = [c for c in coefficients if c.is_cellwise_constant()]
            coefficients = [c for c in coefficients if not c.is_cellwise_constant()]
            integrands = compileUFL(self.integrands, constants=constants, coefficients=coefficients, tempVars=self.tempVars)
        else:
            integrands = self.integrands

        code = [Include('config.h')]
        code += [Include(i) for i in self.gridIncludes]
        #code.append(Include("dune/fem/misc/boundaryidprovider.hh"))

        code += integrands.includes()
        code.append(Include("dune/python/pybind11/pybind11.h"))
        code.append(Include("dune/python/pybind11/extensions.h"))
        code.append(Include("dune/fempy/py/grid/gridpart.hh"))

        if integrands._coefficients:
            code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
        code.append(Include("dune/fempy/py/integrands.hh"))

        nameSpace = NameSpace('Integrands_' + integrands.signature())
        nameSpace.append(integrands.code())
        code.append(nameSpace)

        code.append(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + self.gridType + ' >'))
        if integrands._coefficients:
            coefficients = ['Dune::FemPy::VirtualizedGridFunction< GridPart, ' + c + ' >' for c in integrands._coefficients]
        else:
            coefficients = []
        code.append(TypeAlias('Integrands', nameSpace.name + '::Integrands< ' + ', '.join(['GridPart'] + coefficients) + ' >'))

        writer = SourceWriter()
        writer.emit(code);

        name = self.name()
        writer.openPythonModule(name)
        writer.emit('auto cls = Dune::FemPy::registerIntegrands< Integrands >( module );')
        writer.emit('cls.def( pybind11::init( [] () { return new Integrands(); } ) );')
        writer.closePythonModule(name)

        source = writer.writer.getvalue()
        writer.close()
        return source


def load(grid, integrands, renumbering=None, tempVars=True):
    if isinstance(integrands, Equation):
        integrands = integrands.lhs - integrands.rhs

    source = Source(grid._typeName, grid._includes, integrands, tempVars=tempVars)
    if isinstance(integrands, Form) and renumbering is None:
        coefficients = set(integrands.coefficients())
        renumbering = dict()
        renumbering.update((c, i) for i, c in enumerate(c for c in coefficients if not c.is_cellwise_constant()))
        renumbering.update((c, i) for i, c in enumerate(c for c in coefficients if c.is_cellwise_constant()))

    from dune.generator import builder
    module = builder.load(source.name(), source, "integrands")
    rangeValueTuple, domainValueTuple = source.valueTuples()
    setattr(module.Integrands, "_domainValueType", domainValueTuple)
    setattr(module.Integrands, "_rangeValueType", rangeValueTuple)
    if (renumbering is not None) and not hasattr(module.Integrands, "_renumbering"):
        module.Integrands._setConstant = module.Integrands.__dict__['setConstant']
        module.Integrands._setCoefficient = module.Integrands.__dict__['setCoefficient']
        setattr(module.Integrands, '_renumbering', renumbering)
        setattr(module.Integrands, 'setConstant', setConstant)
        setattr(module.Integrands, 'setCoefficient', setCoefficient)
    return module
