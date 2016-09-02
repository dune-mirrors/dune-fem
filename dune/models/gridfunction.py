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

def gridFunction(grid, code):
    start_time = timeit.default_timer()
    compilePath = os.path.join(os.path.dirname(__file__), '../generated')

    dimRange = grid.dimensionworld
    if not isinstance(grid, types.ModuleType):
        grid = grid._module

    myCodeHash = hashlib.md5(code.encode('utf-8')).hexdigest()
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
            writer.emit('#include <dune/fem/space/common/interpolate.hh>')
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

            writer.openConstMethod('void evaluate', args=['const PointType &x', 'RangeType &value'], targs=['class PointType'])
            writer.emit('const DomainType xGlobal = entity().geometry().global( Dune::Fem::coordinate( x ) );')
            writer.emit(code)
            writer.closeConstMethod()

            writer.openConstMethod('void jacobian', args=['const PointType &x', 'JacobianRangeType &value'], targs=['class PointType'])
            writer.emit('// not implemented')
            writer.closeConstMethod()

            writer.openConstMethod('void hessian', args=['const PointType &x', 'HessianRangeType &value'], targs=['class PointType'])
            writer.emit('// not implemented')
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
            writer.emit('module.def( "get", [] ( const std::string name, const GridPartType &gridPart ) {')
            writer.emit('        LocalFunction local;')
            writer.emit('        return new GridFunction(name, local, gridPart );')
            writer.emit('}, pybind11::keep_alive< 0, 1 >());')
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
