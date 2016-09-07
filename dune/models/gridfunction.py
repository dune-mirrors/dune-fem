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
# import dune.fem.gridpart as gridpart
from dune.fem.gridpart import gridFunctions

# method to add to gridpart.function call
def generatedFunction(grid, name, order, code):
    gf = gridFunction(grid, code)
    return gf.get(name, order, grid)
gridFunctions.update( {"code":generatedFunction} )
def UFLFunction(grid, name, order, expr):
    import ufl
    import dune.models.elliptic as generate
    R = len(expr)
    D = grid.dimension
    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R,), expr), False))
    evaluate = code.replace("result","value")
    jac = []
    for r in range(R):
        jac_form = [\
            ufl.algorithms.expand_indices(ufl.algorithms.expand_derivatives(ufl.algorithms.expand_compounds(\
                ufl.grad(expr)[r,d]*ufl.dx\
            ))) for d in range(D)]
        jac.append( [jac_form[d].integrals()[0].integrand() if not jac_form[d].empty() else 0 for d in range(D)] )
    jac = ufl.as_matrix(jac)
    code = '\n'.join(c for c in generate.generateCode({}, generate.ExprTensor((R,D), jac), False))
    jacobian = code.replace("result","value")
    return generatedFunction(grid,name,order,{"evaluate":evaluate,"jacobian":jacobian})
gridFunctions.update( {"ufl":UFLFunction} )

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

def gridFunction(grid, code):
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

            writer.openStruct(locname, targs=(['class GridPart'] + ['class Range']+ ['class... Coefficients']), bases=(["Dune::Fem::LocalFunctionAdapterHasInitialize"]) )
            writer.typedef(locname + '< GridPart, Range >', 'LocalFunction')
            writer.typedef('GridPart', 'GridPartType')
            writer.typedef('Range', 'RangeType')
            writer.typedef('Dune::Fem::LocalFunctionAdapter< LocalFunction >', 'GridFunction')
            writer.typedef('typename GridPart::template Codim< 0 >::EntityType', 'EntityType')
            writer.typedef('typename EntityType::Geometry::LocalCoordinate', 'LocalCoordinateType')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPart::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename FunctionSpaceType::JacobianRangeType', 'JacobianRangeType')
            writer.typedef('typename FunctionSpaceType::DomainType', 'DomainType')
            writer.typedef('typename FunctionSpaceType::HessianRangeType', 'HessianRangeType')

            writer.openConstMethod('bool init', args=['const EntityType &entity'])
            writer.emit('entity_ = &entity;')
            writer.emit('return true;')
            writer.closeConstMethod()

            writer.openConstMethod('const EntityType &entity')
            writer.emit('return *entity_;')
            writer.closeConstMethod()

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

            writer.section('private')
            writer.openConstMethod('void initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
            writer.emit('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );')
            writer.closeConstMethod()
            writer.emit('')
            writer.emit('mutable const EntityType *entity_ = nullptr;')
            writer.emit('mutable std::tuple< Coefficients... > coefficients_;')
            writer.closeStruct()

            writer.emit('')
            writer.typedef(grid._typeName, 'GridPartType')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPartType::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename FunctionSpaceType::RangeType', 'RangeType')
            writer.typedef(locname + '< GridPartType, RangeType >', 'LocalFunction')
            writer.typedef('Dune::Fem::LocalFunctionAdapter< LocalFunction >', 'GridFunction')

            writer.openNameSpace('Dune')
            writer.openNameSpace('FemPy')
            writer.openPythonModule(pyname)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            writer.emit('pybind11::class_< GridFunction > cls = registerGridFunction< GridFunction >( module, "GridFunction" );')
            writer.emit('module.def( "get", [] ( const std::string name, int order, const GridPartType &gridPart ) {')
            writer.emit('        return new GridFunction(name, LocalFunction(), gridPart, order );')
            writer.emit('}, pybind11::keep_alive< 0, 3 >());')
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
