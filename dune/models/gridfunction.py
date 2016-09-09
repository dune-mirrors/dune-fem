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
def generatedFunction(grid, name, order, code, *args, coefficients=None):
    gf = gridFunction(grid, code, *args, coefficients=coefficients)
    return gf.get(name, order, grid)

gridFunctions.update( {"code" : generatedFunction} )

def UFLFunction(grid, name, order, expr):
    import ufl
    import dune.models.elliptic as generate
    R = len(expr)
    D = grid.dimension
    try:
        _, c = ufl.algorithms.analysis.extract_arguments_and_coefficients(expr)
        coef = set(c)
    except:
        coef = None
    idxConst = 0
    idxCoeff = 0
    for coefficient in coef:
        if coefficient.is_cellwise_constant():
            idx = idxConst
            idxConst += 1
        else:
            idx = idxCoeff
            idxCoeff += 1
        setattr(coefficient, "number", idx)
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
    return generatedFunction(grid, name, order, {"evaluate" : evaluate, "jacobian" : jacobian}, coefficients=coef)

gridFunctions.update( {"ufl" : UFLFunction} )

def dimRangeSplit(code):
    """find the dimRange using @dimrange or counting values
    """
    cpp_code = ''
    codeA = code.split("\n")
    if '@dimrange' in code or '@range' in code:
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

def gridFunction(grid, code, *args, coefficients=None):
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
    wrappername = 'GridFunction_' + myCodeHash + '_' + grid._typeHash

    base = BaseModel(dimRange, myCodeHash)
    if coefficients:
        if type(coefficients) is set:
            for coefficient in coefficients:
                if coefficient.is_cellwise_constant():
                    field = None  # must be improved for 'complex'
                    dimR = 1 if coefficient.ufl_shape==() else coefficient.ufl_shape[0]
                else:
                    field = coefficient.ufl_function_space().ufl_element().field()
                    dimR = coefficient.ufl_shape[0]
                base.coefficients.append({ \
                    'number' : coefficient.number, \
                    'counter' : coefficient.count(), \
                    'dimRange' : dimR,\
                    'constant' : coefficient.is_cellwise_constant(),
                    'field': field } )
    elif args:
        assert len(args) == 1,\
           "arg needs to be a single dict object."
        for arg in args:
            if type(arg) is dict:
                coefNum = 0
                constNum = 0
                for key, value in arg.items():
                    if value[1] == True:
                        num = constNum
                        constNum += 1
                        field = None
                    else:
                        num = coefNum
                        coefNum += 1
                        field = 'double'
                    base.coefficients.append({ \
                        'number' : num, \
                        'dimRange' : value[0], \
                        'constant' : value[1], \
                        'field' : field } )
            else:
                print('Error: arg needs to be dict object containing coefficients.')
                exit(1)

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
            writer.emit('')
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
                writer.typedef(locname + '< GridPart, RangeType, ' + ', '.join(\
                [('Dune::FemPy::VirtualizedLocalFunction< GridPart,'+\
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
            writer.emit(wrappername + '( const std::string name, const GridPart &gridPart, int order ) :')
            writer.emit('    BaseType(name, localFunctionImpl_, gridPart, order) {}')
            writer.emit('LocalFunction& impl() { return localFunctionImpl_; }')
            writer.emit('LocalFunction localFunctionImpl_;')
            writer.closeStruct()
            writer.typedef(wrappername + '< GridPart, RangeType >', 'GFWrapper')

            if base.coefficients:
                base.setCoef(writer, modelClass='LocalFunction', wrapperClass='GFWrapper')

            writer.openNameSpace('Dune')
            writer.openNameSpace('FemPy')
            writer.openPythonModule(pyname)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            writer.emit('pybind11::class_< GFWrapper > cls = registerGridFunction< GFWrapper >( module, "GFWrapper" );')
            if base.coefficients:
                writer.emit('cls.def( "setCoefficient", defSetCoefficient( std::make_index_sequence< std::tuple_size<Coefficients>::value >() ) );')
                writer.emit('cls.def( "setConstant", defSetConstant( std::make_index_sequence< std::tuple_size<typename LocalFunction::ConstantsTupleType>::value >() ) );')
            writer.emit('module.def( "get", [] ( const std::string name, int order, const GridPart &gridPart ) {')
            writer.emit('        return new GFWrapper(name, gridPart, order);')
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
