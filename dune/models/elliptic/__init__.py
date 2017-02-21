from __future__ import absolute_import, division, print_function, unicode_literals

import importlib
import timeit
import types

from ufl.equation import Equation

from dune.source.cplusplus import NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter

from .model import EllipticModel
from .ufl import compileUFL

def generateModel(grid, model, *args, **kwargs):
    start_time = timeit.default_timer()

    if isinstance(model, Equation):
        model = compileUFL(model, *args, **kwargs)

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'ellipticmodel_' + model.signature + "_" + grid._moduleName

    writer = SourceWriter()

    writer.emit(["#include <" + i + ">" for i in grid._includes])
    writer.emit('')
    writer.emit('#include <dune/fem/misc/boundaryidprovider.hh>')
    writer.emit('')
    writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
    writer.emit('#include <dune/corepy/pybind11/extensions.h>')
    writer.emit('')
    writer.emit('#include <dune/fempy/py/grid/gridpart.hh>')
    if model.coefficients:
        writer.emit('#include <dune/fempy/function/virtualizedgridfunction.hh>')
        writer.emit('')
    writer.emit('#include <dune/fem/schemes/diffusionmodel.hh>')

    code = []

    nameSpace = NameSpace("ModelImpl_" + model.signature)
    nameSpace.append(model.code())
    code.append(nameSpace)

    code += [TypeAlias("GridPart", "typename Dune::FemPy::GridPart< " + grid._typeName + " >")]

    rangeTypes = ["Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " >" for c in model.coefficients if not c['constant']]
    coefficients = ["Dune::FemPy::VirtualizedLocalFunction< GridPart, " + r + " >" for r in rangeTypes]
    code += [TypeAlias("Model", nameSpace.name + "::Model< " + ", ".join(["GridPart"] + coefficients) + " >")]

    code += [TypeAlias("ModelWrapper", "DiffusionModelWrapper< Model >"),
             TypeAlias("ModelBase", "typename ModelWrapper::Base")]

    writer.emit(code)

    if model.coefficients:
        model.setCoef(writer)

    writer.openPythonModule(name)
    writer.emit('// export abstract base class')
    writer.emit('if( !pybind11::already_registered< ModelBase >() )')
    writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
    writer.emit('')
    writer.emit('// actual wrapper class for model derived from abstract base')
    writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
    writer.emit('')
    model.export(writer, 'Model', 'ModelWrapper')
    writer.closePythonModule(name)

    if "header" in kwargs:
        with open(kwargs["header"], 'w') as modelFile:
            modelFile.write(writer.writer.getvalue())
    return writer, name


def importModel(grid, model, *args, **kwargs):
    from dune.generator import builder
    if isinstance(model, str):
        with open(model, 'r') as modelFile:
            data = modelFile.read()
        name = data.split('PYBIND11_PLUGIN( ')[1].split(' )')[0]
        builder.load(name, data, "ellipticModel")
        return importlib.import_module("dune.generated." + name)
    writer, name = generateModel(grid, model, *args, **kwargs)
    builder.load(name, writer.writer.getvalue(), "ellipticModel")
    writer.close()
    return importlib.import_module("dune.generated." + name)
