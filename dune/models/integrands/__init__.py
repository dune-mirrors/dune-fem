from __future__ import division, print_function, unicode_literals

from ufl.equation import Equation

from dune.source.cplusplus import Include, NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter

from .model import Integrands
from .ufl import compileUFL

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


def load(grid, integrands, renumbering=None, tempVars=True):
    from dune.common.hashit import hashIt

    if isinstance(integrands, (Form, Equation)):
        integrands, renumbering = compileUFL(integrands, tempVars=tempVars)

    name = 'integrands_' + integrands.signature + '_' + hashIt(grid._typeName)

    code = [Include('config.h')]
    code += [Include(i) for i in grid._includes]
    #code.append(Include("dune/fem/misc/boundaryidprovider.hh"))

    code += integrands.includes()
    code.append(Include("dune/python/pybind11/pybind11.h"))
    code.append(Include("dune/python/pybind11/extensions.h"))
    code.append(Include("dune/fempy/py/grid/gridpart.hh"))

    if integrands._coefficients:
        code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
    code.append(Include("dune/fempy/py/integrands.hh"))

    nameSpace = NameSpace('Integrands_' + integrands.signature)
    nameSpace.append(integrands.code())
    code.append(nameSpace)

    code.append(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + grid._typeName + ' >'))
    if integrands._coefficients:
        rangeTypes = ['Dune::FieldVector< ' + SourceWriter.cpp_fields(c['field']) + ', ' + str(c['dimRange']) + ' >' for c in integrands._coefficients]
        coefficients = ['Dune::FemPy::VirtualizedGridFunction< GridPart, ' + r + ' >' for r in rangeTypes]
    else:
        coefficients = []
    code.append(TypeAlias('Integrands', nameSpace.name + '::Integrands< ' + ', '.join(['GridPart'] + coefficients) + ' >'))

    writer = SourceWriter()
    writer.emit(code);

    writer.openPythonModule(name)
    writer.emit('auto cls = Dune::FemPy::registerIntegrands< Integrands >( module );')
    writer.emit('cls.def( pybind11::init( [] () { return new Integrands(); } ) );')
    writer.closePythonModule(name)

    source = writer.writer.getvalue()
    writer.close()

    from dune.generator import builder
    module = builder.load(name, source, "integrands")
    setattr(module.Integrands, "_domainValueType", integrands.domainValueTuple())
    setattr(module.Integrands, "_rangeValueType", integrands.rangeValueTuple())
    if (renumbering is not None) and not hasattr(module.Integrands, "_renumbering"):
        module.Integrands._setConstant = module.Integrands.__dict__['setConstant']
        module.Integrands._setCoefficient = module.Integrands.__dict__['setCoefficient']
        #setattr(module.Integrands, "_setConstant", getattr(module.Integrands, "setConstant"))
        #setattr(module.Integrands, "_setCoefficient", getattr(module.Integrands, "setCoefficient"))
        setattr(module.Integrands, '_renumbering', renumbering)
        setattr(module.Integrands, 'setConstant', setConstant)
        setattr(module.Integrands, 'setCoefficient', setCoefficient)
    return module


def create(grid, integrands, renumbering=None, tempVars=True):
    return load(grid, integrands, renumbering=renumbering, tempVars=tempVars).Integrands()
