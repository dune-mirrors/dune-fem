from __future__ import absolute_import, division, print_function, unicode_literals

import importlib
import timeit
import types

from ufl import Coefficient, Form
from ufl.equation import Equation

from dune.common.compatibility import isString
from dune.common.hashit import hashIt

from dune.source.cplusplus import NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter

from dune.ufl import GridCoefficient

from .model import EllipticModel
from .ufl import compileUFL


def initModel(model, *args, **kwargs):
    coefficients = kwargs.pop('coefficients', dict())
    if len(args) == 1 and isinstance(args[0], dict):
        coefficients.update(args[0])
        args = []
    else:
        args = list(args)

    coefficientNames = model._coefficientNames
    if len(args) > len(coefficientNames):
        raise ArgumentError('Too many coefficients passed.')
    args += [None] * (len(coefficientNames) - len(args))

    for name, value in kwargs:
        i = coefficientNames.get(name)
        if i is None:
            raise ArgumentError('No such coefficent: ' + name + '.')
        if args[i] is not None:
            raise ArgumentError('Coefficient already given as positional argument: ' + name + '.')
        args[i] = value

    for key, value in coefficients.items():
        if isinstance(key, Coefficient):
            try:
                i = model._renumbering[key]
            except AttributeError:
                raise ArgumentError('Cannot map UFL coefficients, because model was not generated from UFL form.')
            except KeyError:
                raise ArgumentError('No such coefficient: ' + str(key) + '.')
        elif isString(key):
            i = coefficientNames.get(name)
            if i is None:
                raise ArgumentError('No such coefficent: ' + name + '.')
        else:
            raise ArgumentError('Expecting keys of coefficient map to be ufl.Coefficient instances.')
        if args[i] is not None:
            raise ArgumentError('Coefficient already given as positional or keyword argument: ' + str(key) + '.')
        args[i] = value

    if hasattr(model, '_renumbering'):
        for c in (k for k in model._renumbering if isinstance(k, GridCoefficient)):
            i = model._renumbering[c]
            if args[i] is None:
                args[i] = c.gf

    if any(arg is None for arg in args):
        missing = [name for name, i in coefficientNames if args[i] is None]
        raise ArgumentError('Missing coefficients: ' + ', '.join(missing) + '.')

    model._init(*args)


def load(grid, model, *args, **kwargs):
    if isinstance(model, (Equation, Form)):
        model, renumbering = compileUFL(model, *args, **kwargs)
    else:
        renumbering = kwargs.get("renumbering")

    name = 'ellipticmodel_' + model.signature + "_" + hashIt(grid._typeName)

    writer = SourceWriter()

    writer.emit('#include <config.h>')
    writer.emit(["#include <" + i + ">" for i in grid._includes])
    writer.emit('')
    writer.emit('#include <dune/fem/misc/boundaryidprovider.hh>')
    writer.emit('')
    writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
    writer.emit('#include <dune/corepy/pybind11/extensions.h>')
    writer.emit('')
    writer.emit('#include <dune/fempy/py/grid/gridpart.hh>')
    if model.hasCoefficients:
        writer.emit('#include <dune/fempy/function/virtualizedgridfunction.hh>')
        writer.emit('')
    writer.emit('#include <dune/fem/schemes/diffusionmodel.hh>')

    code = []

    nameSpace = NameSpace("ModelImpl_" + model.signature)
    nameSpace.append(model.code())
    code.append(nameSpace)

    code += [TypeAlias("GridPart", "typename Dune::FemPy::GridPart< " + grid._typeName + " >")]

    rangeTypes = ["Dune::FieldVector< " + SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " >" for c in model._coefficients]
    coefficients = ["Dune::FemPy::VirtualizedLocalFunction< GridPart, " + r + " >" for r in rangeTypes]
    code += [TypeAlias("Model", nameSpace.name + "::Model< " + ", ".join(["GridPart"] + coefficients) + " >")]

    code += [TypeAlias("ModelWrapper", "DiffusionModelWrapper< Model >"),
             TypeAlias("ModelBase", "typename ModelWrapper::Base")]

    writer.emit(code)

    if model.hasConstants:
        model.exportSetConstant(writer)

    writer.openPythonModule(name)
    writer.emit('// export abstract base class')
    writer.emit('if( !pybind11::already_registered< ModelBase >() )')
    writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
    writer.emit('')
    writer.emit('// actual wrapper class for model derived from abstract base')
    writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
    writer.emit('')
    for n, number in model._constantNames.items():
        writer.emit('cls.def_property( "' + n + '", ' +
          '[] ( ModelWrapper &self ) { return self.impl().template constant<' + str(number) + '>(); }, ' +
          '[] ( ModelWrapper &self, typename ModelWrapper::Impl::ConstantType<' + str(number) + '>& value) { self.impl().template constant<' + str(number) + '>() = value; }' +
          ');')
    writer.emit('')

    model.export(writer, 'Model', 'ModelWrapper')
    writer.closePythonModule(name)

    source = writer.writer.getvalue()
    writer.close()

    if "header" in kwargs:
        with open(kwargs["header"], 'w') as modelFile:
            modelFile.write(source)

    from dune.generator import builder
    module = builder.load(name, source, "ellipticModel")
    if renumbering is not None:
        setattr(module.Model, '_renumbering', renumbering)
        setattr(module.Model, '_coefficientNames', {c['name']: i for i, c in enumerate(model._coefficients)})
        module.Model._init = module.Model.__dict__['__init__']
        setattr(module.Model, '__init__', initModel)
    return module


def create(grid, model, *args, **kwargs):
    module = load(grid, integrands, renumbering=renumbering, tempVars=tempVars)
    coefficients = kwargs.get('coefficients')
    if coefficients is not None:
        return module.Model(coefficients=coefficients)
    else:
        return module.Model()
