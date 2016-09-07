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
def generatedFunction(grid, name, code, coefficients=None):
    gf = gridFunction(grid, code, coefficients=coefficients)
    return gf.get(name, grid)

gridFunctions.update( {"code" : generatedFunction} )

def UFLFunction(grid, name, expr):
    import ufl
    import dune.models.elliptic as generate
    R = len(expr)
    D = grid.dimension
    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R, ), expr), False))
    evaluate = code.replace("result", "value")
    jac = []
    for r in range(R):
        jac_form = [\
            ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(\
                ufl.grad(expr)[r, d]*ufl.dx\
            ))) for d in range(D)]
        jac.append( [jac_form[d].integrals()[0].integrand() if not jac_form[d].empty() else 0 for d in range(D)] )
    jac = ufl.as_matrix(jac)
    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R, D), jac), False))
    jacobian = code.replace("result", "value")
    try:
        coef = set(expr.coefficients())
    except:
        coef = None
    return generatedFunction(grid, name, {"evaluate" : evaluate, "jacobian" : jacobian}, coefficients=coef)

gridFunctions.update( {"ufl" : UFLFunction} )

def dimRangeSplit(code):
    """find the dimRange using @dimrange or counting values
    """
    cpp_code = ''
    codeA = code.split("\n")
    if '@dimrange' in code or '@range' in code:
        print('@dimrange specified')
        for c in codeA:
            if '@dimrange' in c or '@range' in c:
                dimRange = int( c.split("=",1)[1] )
            else:
                cpp_code += c + "\n"
        cpp_code = cpp_code[:-2]
    else:
        codeB = [c.split("=") for c in codeA]
        codeC = [c[0] for c in codeB if "value" in c[0]]
        dimRange = max( [int(c.split("[")[1].split("]")[0]) for c in codeC] ) + 1
        cpp_code = code
    return ( cpp_code, dimRange )

def gridFunction(grid, code, coefficients=None):
    start_time = timeit.default_timer()
    compilePath = os.path.join(os.path.dirname(__file__), '../generated')

    if type(code) is not dict: code = { 'eval': code }
    cpp_code = ''
    eval = ''
    jac = ''
    hess = ''
    for key, value in code.items():
        if key == 'eval' or key == 'evaluate':
            eval, dimRange = dimRangeSplit(value)
            cpp_code += eval + '\n'
        elif key == 'jac' or key == 'jacobian':
            jac, dimRange = dimRangeSplit(value)
            cpp_code += jac + '\n'
        elif key == 'hess' or key == 'hessian':
            hess, dimRange = dimRangeSplit(value)
            cpp_code += hess + '\n'
        else:
            print(key, ' is not a valid key. Use "eval", "jac" or "hess"')
            exit(1)

    if not isinstance(grid, types.ModuleType):
        grid = grid._module

    myCodeHash = hashlib.md5(cpp_code.encode('utf-8')).hexdigest()
    locname = 'LocalFunction_' + myCodeHash + '_' + grid._typeHash
    pyname = 'localfunction_' + myCodeHash + '_' + grid._typeHash

    base = BaseModel(dimRange, myCodeHash)
    if coefficients:
        base.coefficients = coefficients

    if comm.rank == 0:
        if not os.path.isfile(os.path.join(compilePath, pyname + '.so')):
            writer = SourceWriter(compilePath + '/gridfunction.hh')
            writer.emit(grid._includes)
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
            writer.typedef(grid._typeName, 'GridPartType')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPartType::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename FunctionSpaceType::RangeType', 'RangeType')

            if base.coefficients:
                writer.typedef(locname + '< GridPartType, RangeType, ' + ', '.join(\
                [('Dune::FemPy::VirtualizedLocalFunction< GridPart,'+\
                    'Dune::FieldVector< ' +\
                    SourceWriter.cpp_fields(coefficient['field']) + ', ' +\
                    str(coefficient['dimRange']) + ' > >') \
                    for coefficient in base.coefficients if not coefficient["constant"]])\
                  + ' >', 'LocalFunction')
            else:
                writer.typedef(locname + '< GridPartType, RangeType >', 'LocalFunction')
            writer.typedef('Dune::Fem::LocalFunctionAdapter< LocalFunction >', 'GridFunction')

            if base.coefficients:
                base.setCoef(writer, modelClass='LocalFunction', wrapperClass='GridFunction')

            writer.openNameSpace('Dune')
            writer.openNameSpace('FemPy')
            writer.openPythonModule(pyname)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            writer.emit('pybind11::class_< GridFunction > cls = registerGridFunction< GridFunction >( module, "GridFunction" );')
            if base.coefficients:
                writer.emit('cls.def( "setCoefficient", defSetCoefficient( std::make_index_sequence< std::tuple_size<Coefficients>::value >() ) );')
                writer.emit('cls.def( "setConstant", defSetConstant( std::make_index_sequence< std::tuple_size<typename Model::ConstantsTupleType>::value >() ) );')
            writer.emit('module.def( "get", [] ( const std::string name, const GridPartType &gridPart ) {')
            writer.emit('        LocalFunction local;')
            writer.emit('        return new GridFunction(name, local, gridPart );')
            writer.emit('}')#, pybind11::keep_alive< 0, 2 >());') #error here?
            writer.closePythonModule(pyname)
            writer.closeNameSpace('FemPy')
            writer.closeNameSpace('Dune')

            writer.close()

            cmake = subprocess.Popen(['cmake', '--build', '../../..', '--target', 'gridfunction'], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, 'gridfunction.so'), os.path.join(compilePath, pyname + '.so'))
            print("Compilation took: " , timeit.default_timer()-start_time , "seconds")

        comm.barrier()
        return importlib.import_module('dune.generated.' + pyname)
