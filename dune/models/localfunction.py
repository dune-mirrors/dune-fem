from __future__ import print_function

import importlib
import hashlib
import os
import subprocess
import sys
import timeit
import types
from dune import comm
from dune.source import SourceWriter
from dune.source import BaseModel
# import dune.fem.gridpart as gridpart
from dune.fem.gridpart import gridFunctions

# method to add to gridpart.function call
def generatedFunction(grid, name, order, code, **kwargs):
    coef = kwargs.pop("coefficients", {})
    Gf = gridFunction(grid, code, coef).GFWrapper
    return Gf(name, order, grid, coef)

gridFunctions.update( {"code" : generatedFunction} )

def UFLFunction(grid, name, order, expr, **kwargs):
    import ufl
    import dune.models.elliptic as generate
    from dune.ufl import GridCoefficient
    R = len(expr)
    D = grid.dimension
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
            name = coefficient.str()
        except:
            name = str(coefficient)
        coefficients.append({ \
                    'name' : name, \
                    'number' : idx, \
                    'counter' : coefficient.count(), \
                    'dimRange' : dimR,\
                    'constant' : coefficient.is_cellwise_constant(),
                    'field': field } )

    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R, ), expr), coefficients, False))
    evaluate = code.replace("result", "value")
    jac = []
    for r in range(R):
        jacForm = [\
            ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(\
                ufl.grad(expr)[r, d]*ufl.dx\
            ))) for d in range(D)]
        jac.append( [jacForm[d].integrals()[0].integrand() if not jacForm[d].empty() else 0 for d in range(D)] )
    jac = ufl.as_matrix(jac)
    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R, D), jac), coefficients, False))
    jacobian = code.replace("result", "value")

    code = {"evaluate" : evaluate, "jacobian" : jacobian}
    Gf = gridFunction(grid, code, coefficients).GFWrapper

    coefficients = kwargs.pop("coefficients", {})
    fullCoeff = {c:c.gf for c in coef if isinstance(c, GridCoefficient)}
    fullCoeff.update(coefficients)
    kwargs["coefficients"] = fullCoeff
    return Gf(name, order, grid, **kwargs)

gridFunctions.update( {"ufl" : UFLFunction} )

def gridFunction(grid, code, coefficients):
    startTime = timeit.default_timer()
    compilePath = os.path.join(os.path.dirname(__file__), '../generated')

    if type(code) is not dict: code = { 'eval': code }
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

    if not isinstance(grid, types.ModuleType):
        grid = grid._module

    myCodeHash = hashlib.md5(cppCode.encode('utf-8')).hexdigest()
    locname = 'LocalFunction_' + myCodeHash + '_' + grid._typeHash
    pyname = 'localfunction_' + myCodeHash + '_' + grid._typeHash
    wrappername = 'GridFunction_' + myCodeHash + '_' + grid._typeHash

    base = BaseModel(dimRange, myCodeHash)
    if isinstance(coefficients, dict):
        eval = base.codeCoefficient(eval, coefficients)
        jac = base.codeCoefficient(jac, coefficients)
        hess = base.codeCoefficient(hess, coefficients)
    else:
        base.coefficients.extend(coefficients)

    if comm.rank == 0:
        if not os.path.isfile(os.path.join(compilePath, pyname + '.so')):
            writer = SourceWriter(compilePath + '/localfunction.hh')
            writer.emit("".join(["#include <" + i + ">\n" for i in grid._includes]))
            writer.emit('')
            writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
            writer.emit('#include <dune/corepy/pybind11/extensions.h>')
            writer.emit('')
            writer.emit('#include <dune/fem/space/common/functionspace.hh>')
            writer.emit('#include <dune/fem/function/common/localfunctionadapter.hh>')
            writer.emit('')
            writer.emit('#include <dune/fempy/py/grid/function.hh>')
            writer.emit('')

            base.pre(writer, name=locname, targs=(['class Range']), bases=(["Dune::Fem::LocalFunctionAdapterHasInitialize"]))
            writer.typedef('typename EntityType::Geometry::LocalCoordinate', 'LocalCoordinateType')

            writer.openConstMethod('void evaluate', args=['const PointType &x', 'RangeType &value'], targs=['class PointType'],implemented=eval)
            if eval:
                eval = eval.strip()
                if 'xGlobal' in eval:
                    writer.emit('const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) );')
                writer.emit(eval.split("\n"))
            writer.closeConstMethod()

            writer.openConstMethod('void jacobian', args=['const PointType &x', 'JacobianRangeType &value'], targs=['class PointType'],implemented=jac)
            if jac:
                jac = jac.strip()
                if 'xGlobal' in jac:
                    writer.emit('const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) );')
                writer.emit(jac.split("\n"))
            writer.closeConstMethod()

            writer.openConstMethod('void hessian', args=['const PointType &x', 'HessianRangeType &value'], targs=['class PointType'], implemented=hess)
            if hess:
                hess = hess.strip()
                if 'xGlobal' in hess:
                    writer.emit('const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) );')
                writer.emit(hess.split("\n"))
            writer.closeConstMethod()

            base.post(writer, name=locname, targs=(['class Range']))

            writer.emit('')
            writer.typedef(grid._typeName, 'GridPart')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPart::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename FunctionSpaceType::RangeType', 'RangeType')

            if base.coefficients:
                writer.typedef(locname + '< GridPart, RangeType ' + ''.join(\
                [(', Dune::FemPy::VirtualizedLocalFunction< GridPart,'+\
                    'Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' > >') \
                    for coefficient in base.coefficients if not coefficient["constant"]])\
                  + ' >', 'LocalFunction')
            else:
                writer.typedef(locname + '< GridPart, RangeType >', 'LocalFunction')
            writer.typedef('Dune::Fem::LocalFunctionAdapter< LocalFunction >', 'GridFunction')

            writer.openStruct(wrappername, targs=(['class GridPart'] + ['class Range']), bases=(['GridFunction']))
            writer.typedef('GridFunction', 'BaseType')
            writer.emit(wrappername + '( const std::string name, int order, const GridPart &gridPart ) :')
            writer.emit('    BaseType(name, localFunctionImpl_, gridPart, order) {}')
            writer.emit('LocalFunction& impl() { return localFunctionImpl_; }')
            writer.emit('LocalFunction localFunctionImpl_;')
            writer.closeStruct()
            writer.typedef(wrappername + '< GridPart, RangeType >', 'GFWrapper')

            if base.coefficients:
                base.setCoef(writer, modelClass='LocalFunction', wrapperClass='GFWrapper')

            writer.openPythonModule(pyname)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            writer.emit('pybind11::class_< GFWrapper > cls = Dune::FemPy::registerGridFunction< GFWrapper >( module, "GFWrapper" );')
            writer.emit('')
            base.export(writer, 'LocalFunction', 'GFWrapper', constrArgs = (('name', 'std::string'), ('order', 'int'), ('gridPart', 'GridPart&')), constrKeepAlive='pybind11::keep_alive<0,3>()')
            writer.emit('')
            writer.closePythonModule(pyname)

            writer.close()

            cmake = subprocess.Popen(['cmake', '--build', '../../..', '--target', 'localfunction'], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, 'localfunction.so'), os.path.join(compilePath, pyname + '.so'))
            print("Compilation took: " , timeit.default_timer()-startTime , "seconds")

        comm.barrier()
        return importlib.import_module('dune.generated.' + pyname)
