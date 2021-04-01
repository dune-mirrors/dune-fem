from __future__ import print_function

import importlib
import hashlib
import os
import sys
import timeit
import types

from dune.common.hashit import hashIt
from dune.generator import builder

from dune.source.cplusplus import Include, Method, UnformattedExpression, UnformattedBlock, Variable, Struct, TypeAlias, Constructor, return_
from dune.source.cplusplus import assign
from dune.source.cplusplus import ListWriter, StringWriter, SourceWriter
from dune.source import BaseModel
from dune.source.fem import declareFunctionSpace
import ufl
from ufl import Coefficient, as_vector, replace
from ufl import checks
from ufl.classes import FloatValue, IntValue
from dune.source.cplusplus import maxEdgeLength, UnformattedExpression,\
       NameSpace
from dune.source.algorithm.extractincludes import extractIncludesFromStatements
from dune.ufl import GridFunction, Constant
from dune.ufl.tensors import ExprTensor
from dune.ufl.codegen import uflSignature, TooHighDerivative

from dune.ufl import codegen

class UFLFunctionSource(codegen.ModelClass):
    version = "v1_1"
    def __init__(self, gridType, gridIncludes, expr,
            name,order,
            tempVars=True, virtualize=True,
            predefined=None):
        if len(expr.ufl_shape) == 0:
            expr = as_vector( [ expr ] )
        dimRange = expr.ufl_shape[0]
        predefined = {} if predefined is None else predefined
        codegen.ModelClass.__init__(self,"UFLLocalFunction", [expr],
          virtualize,dimRange=dimRange, predefined=predefined)
        self.evalCode = []
        self.jacCode = []
        self.hessCode = []
        self.gridType = gridType
        self.gridIncludes = gridIncludes

        self.functionName = name
        self.functionOrder = order
        self.expr = expr
        self.tempVars = tempVars
        self.virtualize = virtualize

        self.codeString = self.toString()

    def includes(self):
        incs = set.union(*[extractIncludesFromStatements(stmts) for stmts in (self.evalCode,)])
        ret = [Include(i) for i in incs]
        ret.sort()
        return ret

    def methods(self,code):
        predefined = {}
        x = ufl.SpatialCoordinate(ufl.triangle) # NOTE: to do get right dimension
        predefined[x] = self.spatialCoordinate('x')
        self.predefineCoefficients(predefined, False)
        codegen.generateMethod(code, self.expr,
            'typename FunctionSpaceType::RangeType', 'evaluate',
            returnResult=False,
            args=['const Point &x'],
            targs=['class Point'], const=True,
            predefined=predefined)
        if checks.is_globally_constant(self.expr):
            code.append( Method('void', 'jacobian', targs=['class Point'],
                args=['const Point &x','typename FunctionSpaceType::JacobianRangeType &result'],
                code=['result=typename FunctionSpaceType::JacobianRangeType(0);'], const=True))
            code.append( Method('void', 'hessian', targs=['class Point'],
                args=['const Point &x','typename FunctionSpaceType::HessianRangeType &result'],
                code=['result=typename FunctionSpaceType::HessianRangeType(0);'], const=True))
        else:
            try:
                codegen.generateMethod(code, ufl.grad(self.expr),
                    'typename FunctionSpaceType::JacobianRangeType', 'jacobian',
                    returnResult=False,
                    args=['const Point &x'],
                    targs=['class Point'], const=True,
                    predefined=predefined)
            except Exception as e:
                code.append( Method('void', 'jacobian', targs=['class Point'],
                    args=['const Point &x','typename FunctionSpaceType::JacobianRangeType &result'],
                    code=['DUNE_THROW(Dune::NotImplemented,"jacobian method could not be generated for local function ('+repr(e)+').");',
                          'result=typename FunctionSpaceType::JacobianRangeType(0);'], const=True))
                pass
            try:
                codegen.generateMethod(code, ufl.grad(ufl.grad(self.expr)),
                    'typename FunctionSpaceType::HessianRangeType', 'hessian',
                    returnResult=False,
                    args=['const Point &x'],
                    targs=['class Point'], const=True,
                    predefined=predefined)
            except Exception as e:
                code.append( Method('void', 'hessian', targs=['class Point'],
                    args=['const Point &x','typename FunctionSpaceType::HessianRangeType &result'],
                    code=['DUNE_THROW(Dune::NotImplemented,"hessian method could not be generated for local function ('+repr(e)+')");',
                          'result=typename FunctionSpaceType::HessianRangeType(0);'], const=True))
                pass

    def _signature(self):
        from dune.common.hashit import hashIt
        coeffTypes = ','.join(self.coefficientCppTypes)
        return hashIt(self.codeString + coeffTypes)
    def signature(self):
        return self._signature()
        return uflSignature(None, *self.coefficientCppTypes,
                            *self._constantNames, self.expr)+UFLFunctionSource.version

    def name(self):
        from dune.common.hashit import hashIt
        if self.virtualize:
            return 'localfunction_' + self._signature() + '_' + hashIt(self.gridType)
        else:
            return 'localfunction_nv' + self._signature() + '_' + hashIt(self.gridType)

    def compileUFL(self):
        pass

    # actual code generation (the generator converts this object to a string)
    def __str__(self):
        coefficients = self.coefficientCppTypes
        self.compileUFL()
        code = [Include('config.h')]
        code += [Include(i) for i in self.gridIncludes]
        code += [Include(i) for i in self.includeFiles]

        code += self.includes()
        code.append(Include("dune/python/pybind11/pybind11.h"))
        code.append(Include("dune/python/pybind11/extensions.h"))
        code.append(Include("dune/fempy/py/grid/gridpart.hh"))
        code.append(Include('dune/common/exceptions.hh'))

        if self._coefficients:
            if self.virtualize:
                code.append(Include("dune/fempy/function/virtualizedgridfunction.hh"))
            else:
                for c in self._coefficients:
                    for i in c._includes:
                        code.append(Include(i))
        code.append(Include("dune/fempy/py/ufllocalfunction.hh"))

        # main part
        nameSpace = NameSpace('UFLLocalFunctions_' + self.signature())
        nameSpace.append(self.codeString)
        code.append(nameSpace)
        gridPartName = 'typename Dune::FemPy::GridPart< ' + self.gridType + ' >'
        localFunctionName = nameSpace.name+'::UFLLocalFunction< ' + ', '.join([gridPartName] + self.coefficientCppTypes) + ' >'

        name = self.name()
        writer = SourceWriter()
        writer.emit(code)
        writer.openPythonModule(name)
        # writer.emit(TypeAlias('GridPart', 'typename Dune::FemPy::GridPart< ' + self.gridType + ' >'))
        # writer.emit(TypeAlias('GridPart', gridPartName))
        writer.emit(TypeAlias('LocalFunctionType', localFunctionName))
        writer.emit('auto cls = Dune::Python::insertClass<LocalFunctionType>(module,"UFLLocalFunction",Dune::Python::GenerateTypeName("'+localFunctionName+'"), Dune::Python::IncludeFiles({"python/dune/generated/'+name+'.cc"})).first;')
        writer.emit('Dune::FemPy::registerUFLLocalFunction( module, cls );')

        initArgs = 'pybind11::object gridView, const std::string &name, int order'
        keepAlive = 'pybind11::keep_alive< 1, 2 >()'
        initCall = 'Dune::FemPy::gridPart<'+self.gridType+'>(gridView),name,order'
        if coefficients:
            coefficientNames = self.coefficientNames
            initArgs += ', ' + ', '.join('const ' + t + ' &' + n for t, n in zip(coefficients, coefficientNames))
            keepAlive += ', ' + ', '.join('pybind11::keep_alive< 1, ' + str(i+3) + ' >()' for i in range(len(coefficientNames)))
            initCall += ', ' + ', '.join(coefficientNames)
        writer.emit('cls.def( pybind11::init( [] ( ' + initArgs + ' ) {'
                + 'return new '+localFunctionName+' ( ' + initCall
                + '); } ), ' + keepAlive + ' );')

        for t, n, sn in zip(self.constantTypes, self.constantNames, self.constantShortNames):
            te = localFunctionName+"::" + t
            writer.emit('cls.def_property( "' + sn + '", [] ( '+localFunctionName+' &self ) -> ' + te + ' { return self.' + n + '(); }, [] ( '+localFunctionName+' &self, const ' + te + ' &v ) { self.' + n + '() = v; } );')
        writer.emit('cls.def_property_readonly( "virtualized", [] ( '+localFunctionName+'& ) -> bool { return '+str(self.virtualize).lower()+';});')

        writer.closePythonModule(name)
        source = writer.writer.getvalue()
        writer.close()

        source = "#ifndef GUARD_"+self.signature()+\
                 "\n#define GUARD_"+self.signature()+"\n"+\
                 source+\
                 "\n#endif\n"

        return source
    def toString(self):
        writer = SourceWriter()
        writer.emit(self.code())
        source = writer.writer.getvalue()
        writer.close()
        return source


def init(lf, gridView, name, order, *args, **kwargs):
    coefficients = kwargs.pop('coefficients', dict())
    coefficientNames = lf._coefficientNames
    if len(args) == 1 and isinstance(args[0], dict):
        coefficients.update(args[0])
    args = []

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
                i = lf._renumbering[key]
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

    if hasattr(lf, '_renumbering'):
        for c, i in lf._renumbering.items():
            if isinstance(c, GridFunction):
                if args[i] is None:
                    args[i] = c.gf

    if any(arg is None for arg in args):
        missing = [name for name, i in coefficientNames.items() if args[i] is None]
        raise ValueError('Missing coefficients: ' + ', '.join(missing) + '.')

    lf.base.__init__(lf,gridView,name,order,*args)
    for c in lf._constants:
        c.registerModel(lf)

def setConstant(lf, index, value):
    try:
        index = lf._renumbering[index]
    except KeyError:
        pass
    lf._setConstant(index, value)

def UFLFunction(grid, name, order, expr, renumbering=None, virtualize=True, tempVars=True,
                predefined=None, **kwargs):
    scalar = False
    if type(expr) == list or type(expr) == tuple:
        expr = ufl.as_vector(expr)
    elif type(expr) == int or type(expr) == float:
        expr = ufl.as_vector( [expr] )
        scalar = True
    try:
        if expr.ufl_shape == ():
            expr = ufl.as_vector([expr])
            scalar = True
    except:
        return None
    _, coeff_ = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
    coeff   = {c : c.toVectorCoefficient()[0] for c in coeff_ if len(c.ufl_shape) == 0 and not c.is_cellwise_constant()}
    expr = replace(expr,coeff)

    if len(expr.ufl_shape) > 1:
        raise AttributeError("can only generate grid functions from vector values UFL expressions not from expressions with shape=",expr.ufl_shape)

    # set up the source class
    source = UFLFunctionSource(grid._typeName, grid._includes, expr,
            name,order,
            tempVars=tempVars,virtualize=virtualize,
            predefined=predefined)

    coefficients = source.coefficientList
    numCoefficients = len(coefficients)
    if renumbering is None:
        renumbering = dict()
        renumbering.update((c, i) for i, c in enumerate(sorted((c for c in coefficients if not c.is_cellwise_constant()), key=lambda c: c.count())))
        renumbering.update((c, i) for i, c in enumerate(c for c in coefficients if c.is_cellwise_constant()))
    coefficientNames = ['coefficient' + str(i) if n is None else n for i, n in enumerate(getattr(c, 'name', None) for c in coefficients if not c.is_cellwise_constant())]

    # call code generator
    from dune.generator import builder
    module = builder.load(source.name(), source, "UFLLocalFunction")

    class LocalFunction(module.UFLLocalFunction):
        def __init__(self, gridView, name, order, *args, **kwargs):
            self.base = module.UFLLocalFunction
            self._coefficientNames = {n: i for i, n in enumerate(source.coefficientNames)}
            if renumbering is not None:
                self._renumbering = renumbering
                self._setConstant = self.setConstant # module.UFLLocalFunction.__dict__['setConstant']
                self.setConstant = lambda *args: setConstant(self,*args)
            self.constantShape = source._constantShapes
            self._constants = [c for c in source.constantList if isinstance(c,Constant)]
            self.scalar = scalar
            init(self, gridView, name, order, *args, **kwargs)

    return LocalFunction
