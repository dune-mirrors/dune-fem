from __future__ import absolute_import, division, print_function, unicode_literals

import importlib
import timeit
import types

from ufl import Coefficient, Form
from ufl.equation import Equation

from dune.common.utility import isString
from dune.common.hashit import hashIt

from dune.source.cplusplus import NameSpace, TypeAlias
from dune.source.cplusplus import SourceWriter

from dune.ufl import GridFunction

from .ufl import compileUFL


def initModel(model, *args, **kwargs):
    coefficients = kwargs.pop('coefficients', dict())
    if len(args) == 1 and isinstance(args[0], dict):
        coefficients.update(args[0])
        args = []
    else:
        args = list(args)

    try:
        coefficientNames = model._coefficientNames
    except:
        coefficientNames = []
    if len(args) > len(coefficientNames):
        raise ValueError('Too many coefficients passed.')
    args += [None] * (len(coefficientNames) - len(args))

    for name, value in kwargs:
        i = coefficientNames.get(name)
        if i is None:
            raise ValueError('No such coefficent: ' + name + '.')
        if args[i] is not None:
            raise ValueError('Coefficient already given as positional argument: ' + name + '.')
        args[i] = value

    for key, value in coefficients.items():
        if isinstance(key, Coefficient):
            try:
                i = model._renumbering[key]
            except AttributeError:
                raise ValueError('Cannot map UFL coefficients, because model was not generated from UFL form.')
            except KeyError:
                raise ValueError('No such coefficient: ' + str(key) + '.')
        elif isString(key):
            i = coefficientNames.get(name)
            if i is None:
                raise ValueError('No such coefficent: ' + name + '.')
        else:
            raise ValueError('Expecting keys of coefficient map to be ufl.Coefficient instances.')
        if args[i] is not None:
            raise ValueError('Coefficient already given as positional or keyword argument: ' + str(key) + '.')
        args[i] = value

    if hasattr(model, '_renumbering'):
        for c in (k for k in model._renumbering if isinstance(k, GridFunction)):
            i = model._renumbering[c]
            if args[i] is None:
                args[i] = c.gf

    if any(arg is None for arg in args):
        missing = [name for name, i in coefficientNames.items() if args[i] is None]
        raise ValueError('Missing coefficients: ' + ', '.join(missing) + '.')

    model._init(*args)


def load(grid, model, *args, modelPatch=[None,None], virtualize=True, **kwargs):
    if not isinstance(modelPatch,list) and not isinstance(modelPatch,tuple):
        modelPatch = [modelPatch,None]

    from dune.generator import builder
    if isinstance(model, (Equation, Form)):
        model = compileUFL(model, modelPatch[1], *args, **kwargs)
        renumbering = model.coefficients.copy()
        renumbering.update(model.constants)
    else:
        renumbering = kwargs.get("renumbering")

    if isinstance(model, str):
        with open(model, 'r') as modelFile:
             data = modelFile.read()
        name = data.split('PYBIND11_MODULE( ')[1].split(',')[0]
        endPos = name.find('_')
        modelName = name[0:endPos]
        module = builder.load(name, data, modelName)
        renumbering = {}
        if renumbering is not None:
            setattr(module.Model, '_renumbering', renumbering)
            module.Model._init = module.Model.__dict__['__init__']
            setattr(module.Model, '__init__', initModel)
        return module

    if modelPatch[0] is not None:
        modelPatch[0](model)
    else:
        modelPatch = None

    signature = ("" if virtualize else "nv") + model.signature + "_" + hashIt(grid.cppTypeName)
    name = model.baseName + '_' + signature

    writer = SourceWriter()

    writer.emit("#ifndef GuardModelImpl_" + signature)
    writer.emit("#define GuardModelImpl_" + signature)
    writer.emit("#define USING_DUNE_PYTHON 1")

    writer.emit('#include <config.h>')
    writer.emit(["#include <" + i + ">" for i in grid.cppIncludes])
    writer.emit('')
    writer.emit('#include <dune/fem/misc/boundaryidprovider.hh>')
    writer.emit('')
    writer.emit('#include <dune/python/pybind11/pybind11.h>')
    writer.emit('#include <dune/python/pybind11/extensions.h>')
    writer.emit('')
    writer.emit('#include <dune/fempy/py/grid/gridpart.hh>')
    if model.hasCoefficients:
        writer.emit('#include <dune/fempy/function/virtualizedgridfunction.hh>')
        writer.emit('')
    if 'virtualModel' in kwargs:
        virtualModel = kwargs.pop('virtualModel')
    else:
        virtualModel = 'dune/fem/schemes/conservationlawmodel.hh'
    writer.emit('#include <' + virtualModel + '>')

    nameSpace = NameSpace("ModelImpl_" + signature)
    if modelPatch:
        nameSpace.append(model.code(model))
    else:
        nameSpace.append(model.code())

    writer.emit(nameSpace)

    writer.openNameSpace("ModelImpl_" + signature)
    gridPartType = "typename Dune::FemPy::GridPart< " + grid.cppTypeName + " >"
    rangeTypes = ["Dune::FieldVector< " +
            SourceWriter.cpp_fields(c['field']) + ", " + str(c['dimRange']) + " >"\
            for c in model._coefficients]
    coefficients = ["Dune::FemPy::VirtualizedGridFunction<"+gridPartType+", " + r + " >"
                    if not c['typeName'].startswith("Dune::Python::SimpleGridFunction") \
                    else c['typeName'] \
            for r,c in zip(rangeTypes,model._coefficients)]
    modelType = nameSpace.name + "::Model< " + ", ".join([gridPartType] + coefficients) + " >"
    if not virtualize:
        wrapperType = modelType
    else:
        wrapperType = model.modelWrapper.replace(" Model ",modelType)
    if model.hasConstants:
        model.exportSetConstant(writer, modelClass=modelType, wrapperClass=wrapperType)

    writer.closeNameSpace("ModelImpl_" + signature)
    writer.openPythonModule(name)
    code = []
    code += [TypeAlias("GridPart", gridPartType)]
    code += [TypeAlias("Model", nameSpace.name + "::Model< " + ", ".join(["GridPart"] + coefficients) + " >")]
    if virtualize:
        modelType = model.modelWrapper
        code += [TypeAlias("ModelWrapper", model.modelWrapper)]
        code += [TypeAlias("ModelBase", "typename ModelWrapper::Base")]
    else:
        modelType = nameSpace.name + "::Model< " + ", ".join([gridPartType] + coefficients) + " >"
        code += [TypeAlias("ModelWrapper", "Model")]
    writer.emit(code)

    if virtualize:
        writer.emit('// export abstract base class')
        writer.emit('if( !pybind11::already_registered< ModelBase >() )')
        writer.emit('  pybind11::class_< ModelBase >( module, "ModelBase" );')
        writer.emit('')
        writer.emit('// actual wrapper class for model derived from abstract base')
        # writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model", pybind11::base< ModelBase >() );')
        writer.emit('auto cls = Dune::Python::insertClass<ModelWrapper,ModelBase>(module,"Model",'+\
                        'Dune::Python::GenerateTypeName("'+modelType+'"),'+\
                        'Dune::Python::IncludeFiles({})).first;')
    else:
        # writer.emit('pybind11::class_< ModelWrapper > cls( module, "Model" );')
        writer.emit('auto cls = Dune::Python::insertClass<ModelWrapper>(module,"Model",'+\
                        'Dune::Python::GenerateTypeName("'+modelType+'"),'+\
                        'Dune::Python::IncludeFiles({"python/dune/generated/'+name+'.cc"})).first;')
    writer.emit('cls.def_property_readonly( "dimRange", [] ( ModelWrapper & ) { return ' + str(model.dimRange) + '; } );')
    hasDirichletBC = 'true' if model.hasDirichletBoundary else 'false'
    writer.emit('cls.def_property_readonly( "hasDirichletBoundary", [] ( ModelWrapper& ) -> bool { return '+hasDirichletBC+';});')
    writer.emit('')
    for n, number in model._constantNames.items():
        writer.emit('cls.def_property( "' + n + '", ' +
          '[] ( ModelWrapper &self ) { return self.template constant<' + str(number) + '>(); }, ' +
          '[] ( ModelWrapper &self, typename ModelWrapper::ConstantType<' + str(number) + '>& value) { self.template constant<' + str(number) + '>() = value; }' +
          ');')
    writer.emit('')

    model.export(writer, 'Model', 'ModelWrapper',nameSpace="ModelImpl_"+signature)
    writer.closePythonModule(name)
    writer.emit("#endif // GuardModelImpl_" + signature)

    source = writer.writer.getvalue()
    writer.close()

    if "header" in kwargs:
        with open(kwargs["header"], 'w') as modelFile:
            modelFile.write(source)

    endPos = name.find('_')
    modelName = name[0:endPos]
    module = builder.load(name, source, modelName)
    if (renumbering is not None) and (module.Model.__dict__['__init__'] != initModel):
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
