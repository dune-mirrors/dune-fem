from __future__ import print_function

import importlib
import hashlib
import os
import subprocess
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

# method to add to gridpart.function call
def generatedFunction(grid, name, order, code, **kwargs):
    coef = kwargs.pop("coefficients", {})
    const = kwargs.pop("constants", {})
    Gf = gridFunction(grid, code, coef, const).GFWrapper
    return Gf(name, order, grid, coef)


def generateCode(predefined, tensor, coefficients, tempVars=True):
    from dune.ufl import codegen
    keys = tensor.keys()
    expressions = [tensor[i] for i in keys]
    preamble, results = codegen.generateCode(predefined, expressions, coefficients, tempVars=tempVars)
    result = Variable('auto', 'result')
    return preamble + [assign(result[i], r) for i, r in zip(keys, results)]


def UFLFunction(grid, name, order, expr, **kwargs):
    import ufl
    from dune.ufl import GridCoefficient
    from dune.ufl.tensors import ExprTensor
    try:
        _, c = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
        coef = set(c)
    except:
        coef = None
    idxConst = 0
    idxCoeff = 0
    coefficients = []
    for coefficient in coef:
        if coefficient.is_cellwise_constant():
            idx = idxConst
            idxConst += 1
            field = None  # must be improved for 'complex'
            dimR = 1 if coefficient.ufl_shape==() else coefficient.ufl_shape[0]
        else:
            idx = idxCoeff
            idxCoeff += 1
            field = coefficient.ufl_function_space().ufl_element().field()
            dimR = coefficient.ufl_shape[0]

        try:
            name = getattr(coefficient, "name")
        except AttributeError:
            name = str(coefficient)
        coefficients.append({ \
                    'name' : name, \
                    'number' : idx, \
                    'counter' : coefficient.count(), \
                    'dimRange' : dimR,\
                    'constant' : coefficient.is_cellwise_constant(),
                    'field': field } )

    if expr.ufl_shape == ():
        expr = ufl.as_vector([expr])
    dimR = expr.ufl_shape[0]
    dimD = grid.dimension

    writer = SourceWriter(ListWriter())
    writer.emit(generateCode({}, ExprTensor((dimR, ), expr), coefficients, False), context=Method('void', 'evaluate'))
    code = '\n'.join(writer.writer.lines)
    evaluate = code.replace("result", "value")
    jac = []
    for r in range(dimR):
        jacForm = [\
            ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(\
                ufl.grad(expr)[r, d]*ufl.dx\
            ))) for d in range(dimD)]
        jac.append( [jacForm[d].integrals()[0].integrand() if not jacForm[d].empty() else 0 for d in range(dimD)] )
    jac = ufl.as_matrix(jac)
    writer = SourceWriter(ListWriter())
    writer.emit(generateCode({}, ExprTensor((dimR, dimD), jac), coefficients, False), context=Method('void', 'jacobian'))
    code = '\n'.join(writer.writer.lines)
    jacobian = code.replace("result", "value")

    code = {"evaluate" : evaluate, "jacobian" : jacobian}
    Gf = gridFunction(grid, code, coefficients, None).GFWrapper

    coefficients = kwargs.pop("coefficients", {})
    fullCoeff = {c:c.gf for c in coef if isinstance(c, GridCoefficient)}
    fullCoeff.update(coefficients)
    kwargs["coefficients"] = fullCoeff
    return Gf(name, order, grid, **kwargs)


def gridFunction(grid, code, coefficients, constants):
    startTime = timeit.default_timer()

    if type(code) is not dict:
        code = {'eval': code}
    cppCode = ''
    eval = ''
    jac = ''
    hess = ''
    for key, value in code.items():
        if key == 'eval' or key == 'evaluate':
            eval, dimRange = BaseModel.codeDimRange(value)
            cppCode += eval + '\n'
        elif key == 'jac' or key == 'jacobian':
            jac, dimRange = BaseModel.codeDimRange(value)
            cppCode += jac + '\n'
        elif key == 'hess' or key == 'hessian':
            hess, dimRange = BaseModel.codeDimRange(value)
            cppCode += hess + '\n'
        else:
            print(key, ' is not a valid key. Use "eval", "jac" or "hess"')
            exit(1)

    # if not isinstance(grid, types.ModuleType):
    #     grid = grid._module

    if isinstance(coefficients, dict):
        for entry in coefficients.items():
            cppCode += entry[0]
    else:
        for coefficient in coefficients:
            cppCode += str(coefficient.get('name'))

    myCodeHash = hashIt(cppCode)
    locname = 'LocalFunction_' + myCodeHash + '_' + hashIt(grid._typeName)
    pyname = 'localfunction_' + myCodeHash + '_' + hashIt(grid._typeName)
    wrappername = 'GridFunction_' + myCodeHash + '_' + hashIt(grid._typeName)

    base = BaseModel(dimRange, myCodeHash)
    if isinstance(coefficients, dict):
        eval = base.codeCoefficient(eval, coefficients, constants)

        jac = base.codeCoefficient(jac, coefficients, constants)
        hess = base.codeCoefficient(hess, coefficients, constants)
    else:
        base.coefficients.extend(coefficients)

    code = []

    code = [Include(i) for i in grid._includes]

    code.append(Include("dune/corepy/pybind11/pybind11.h"))
    code.append(Include("dune/corepy/pybind11/extensions.h"))
    code.append(Include("dune/fem/space/common/functionspace.hh"))
    code.append(Include("dune/fem/function/common/localfunctionadapter.hh"))

    code.append(Include("dune/fempy/py/grid/gridpart.hh"))
    code.append(Include("dune/fempy/py/grid/function.hh"))

    struct = Struct(locname, targs=['class GridPart', 'class Range', 'class... Coefficients'])
    struct.append(base.pre(name=locname))

    struct.append(TypeAlias('LocalCoordinateType', 'typename EntityType::Geometry::LocalCoordinate'))

    evaluate = Method('void', 'evaluate', args=['const Point &x', 'RangeType &value'], targs=['class Point'], const=True)
    if eval:
        eval = eval.strip()
        if 'xGlobal' in eval:
            evaluate.append(UnformattedExpression('void', 'const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) )'))
        evaluate.append(UnformattedBlock(eval.split('\n')))
    else:
        evaluate.append(UnformattedBlock('DUNE_THROW( Dune::NotImplemented, "evaluate not implemented." );'));
    struct.append(evaluate)

    jacobian = Method('void', 'jacobian', args=['const Point &x', 'JacobianRangeType &value'], targs=['class Point'], const=True)
    if jac:
        jac = jac.strip()
        if 'xGlobal' in jac:
            jacobian.append(UnformattedExpression('void', 'const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) )'))
        jacobian.append(UnformattedBlock(jac.split('\n')))
    else:
        jacobian.append(UnformattedBlock('DUNE_THROW( Dune::NotImplemented, "jacobian not implemented." );'));
    struct.append(jacobian)

    hessian = Method('void', 'hessian', args=['const Point &x', 'HessianRangeType &value'], targs=['class Point'], const=True)
    if hess:
        hess = hess.strip()
        if 'xGlobal' in hess:
            hessian.append(UnformattedExpression('void', 'const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) )'))
        hessian.append(UnformattedBlock(hess.split('\n')))
    else:
        hessian.append(UnformattedBlock('DUNE_THROW( Dune::NotImplemented, "hessian not implemented." );'));
    struct.append(hessian)

    struct.append(base.post())
    code.append(struct)

    code.append(TypeAlias('GridPartType', 'Dune::FemPy::GridPart< ' + grid._typeName + ' >'))
    code.append(TypeAlias('GridView', 'typename GridPartType::GridViewType'))
    code.append(declareFunctionSpace("typename GridPartType::ctype", "double", UnformattedExpression("int", "GridPartType::dimensionworld"), dimRange))

    if base.coefficients:
        code.append(TypeAlias('LocalFunction', locname + '< GridPartType, RangeType '
            + ''.join([(', Dune::FemPy::VirtualizedLocalFunction< GridPartType, Dune::FieldVector< ' \
                      + SourceWriter.cpp_fields(coefficient['field']) + ', ' + str(coefficient['dimRange']) + ' > >') \
                      for coefficient in base.coefficients if not coefficient["constant"]]) + ' >'))
    else:
        code.append(TypeAlias('LocalFunction', locname + '< GridPartType, RangeType >'))
    code.append(TypeAlias('GridFunction', 'Dune::Fem::LocalFunctionAdapter< LocalFunction >'))

    wrapper = Struct(wrappername, targs=['class GridView', 'class Range'], bases=['GridFunction'])
    wrapper.append(TypeAlias('BaseType', 'GridFunction'))
    wrapper.append(Constructor(args=['const std::string name', 'int order', 'pybind11::handle gridView'], init=['BaseType(name, LocalFunction(), Dune::FemPy::gridPart< GridView >( gridView ), order )']))
    wrapper.append(Method('LocalFunction &', 'impl', code=return_(UnformattedExpression('auto', 'this->localFunctionImpl()'))))
    code.append(wrapper)

    code.append(TypeAlias('GFWrapper', wrappername + '< GridView, RangeType >'))

    writer = SourceWriter(StringWriter())
    writer.emit(code)

    if base.coefficients:
        base.setCoef(writer, modelClass='LocalFunction', wrapperClass='GFWrapper')

    writer.openPythonModule(pyname)
    writer.emit('')
    writer.emit('// export function class')
    writer.emit('')
    writer.emit('pybind11::class_< GFWrapper > cls = Dune::FemPy::registerGridFunction< GFWrapper >( module, "GFWrapper" );')
    writer.emit('')
    base.export(writer, 'LocalFunction', 'GFWrapper', constrArgs = (('name', 'std::string'), ('order', 'int'), ('gridView', 'pybind11::handle')), constrKeepAlive='pybind11::keep_alive<0,3>()' )
    writer.emit('')
    writer.closePythonModule(pyname)

    builder.load(pyname, writer.writer.getvalue(), "localFunction")
    writer.close()

    return importlib.import_module('dune.generated.' + pyname)
