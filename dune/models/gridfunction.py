from __future__ import print_function

import importlib
import os
import subprocess
import sys
import types
from dune import comm
from dune.source import SourceWriter

def gridFunction(grid, code):
    compilePath = os.path.join(os.path.dirname(__file__), '../generated')

    dimRange = grid.dimensionworld
    if not isinstance(grid, types.ModuleType):
        grid = grid._module

    gfnumber = 0
    check = 0
    while check == 0:
        name = 'gridfunction_' + str(gfnumber) + '_' + grid._typeHash
        if os.path.isfile(os.path.join(compilePath, name + '.so')):
            gfnumber += 1
        else:
            check = 1
    gfname = 'GridFunction_' + str(gfnumber) + '_' + grid._typeHash
    locname = 'LocalFunction_' + str(gfnumber) + '_' + grid._typeHash

    if comm.rank == 0:
        if not os.path.isfile(os.path.join(compilePath, name + '.so')):
            writer = SourceWriter(compilePath + '/gridfunction.hh')
            writer.emit(grid._includes)
            writer.emit('')
            writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
            writer.emit('#include <dune/corepy/pybind11/extensions.h>')
            writer.emit('')
            writer.emit('#include <dune/fem/space/common/functionspace.hh>')
            writer.emit('#include <dune/fem/space/common/interpolate.hh>')
            writer.emit('#include <dune/fempy/py/grid/function.hh>')
            writer.emit('')

            writer.openStruct(locname, targs=(['class GridPart'] + ['class Range']))
            writer.typedef('GridPart', 'GridPartType')
            writer.typedef('Range', 'RangeType')
            writer.typedef('typename GridPart::template Codim< 0 >::EntityType', 'EntityType')
            writer.typedef('typename EntityType::Geometry::LocalCoordinate', 'LocalCoordinateType')

            writer.openConstMethod('const EntityType &entity')
            writer.emit('return *entity_;')
            writer.closeConstMethod()

            writer.openConstMethod('void evaluate', args=['const PointType &x', 'RangeType &value'], targs=['class PointType'])
            writer.emit(code)
            writer.closeConstMethod()

            writer.openConstMethod('void jacobian' args=['const PointType &x', 'RangeType &value'], targs=['class PointType'])
            writer.emit('// not implemented')
            writer.closeConstMethod

            writer.emit('mutable const EntityType *entity_ = nullptr;')
            writer.closeStruct()


            writer.openStruct(gfname, targs=(['class GridPart'] + ['class Range'] + ['class... Coefficients']))

            writer.typedef('GridPart', 'GridPartType')
            writer.typedef('Range', 'RangeType')
            writer.typedef(locname + '< GridPart, Range >', 'LocalFunctionType')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPart::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename GridPart::template Codim< 0 >::EntityType', 'EntityType')
            writer.typedef('typename FunctionSpaceType::DomainType', 'DomainType')
            writer.typedef('typename EntityType::Geometry::LocalCoordinate', 'LocalCoordinateType')
            writer.typedef('typename FunctionSpaceType::JacobianRangeType', 'JacobianRangeType')
            writer.typedef('typename FunctionSpaceType::HessianRangeType', 'HessianRangeType')

            writer.openConstMethod('bool init', args=['const EntityType &entity'] + ['const GridPartType &gridPart'] + ['const LocalFunctionType &localFunction'])
            writer.emit('entity_ = &entity;')
            writer.emit('gridPart_ = &gridPart;')
            writer.emit('localFunction_ = &localFunction;')
            writer.emit('return true;')
            writer.closeConstMethod()

            writer.openConstMethod('const EntityType &entity')
            writer.emit('return *entity_;')
            writer.closeConstMethod()

            writer.openConstMethod('GridPartType &gridPart')
            writer.emit('return *gridPart_;')
            writer.closeConstMethod()

            writer.openConstMethod('LocalFunctionType &localFunction', args=['const EntityType &entity'])
            writer.emit('return *localFunction_;')
            writer.closeConstMethod()

            writer.openConstMethod('std::string name')
            writer.emit('return "' + name + '";')
            writer.closeConstMethod()

            writer.section('private')
            writer.openConstMethod('void initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
            writer.emit('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );')
            writer.closeConstMethod()
            writer.emit('')
            writer.emit('mutable const EntityType *entity_ = nullptr;')
            writer.emit('mutable const GridPartType *gridPart_ = nullptr;')
            writer.emit('mutable const LocalFunctionType *localFunction_ = nullptr;')
            writer.emit('mutable std::tuple< Coefficients... > coefficients_;')
            writer.closeStruct()

            writer.emit('')
            writer.typedef(grid._typeName, 'GridPart')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPart::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename FunctionSpaceType::RangeType', 'RangeType')
            writer.typedef(gfname + '< GridPart, RangeType >', 'GridFunction')

            writer.openNameSpace('Dune')
            writer.openNameSpace('FemPy')
            writer.openPythonModule(name)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            writer.emit('pybind11::class_< GridFunction > cls = registerGridFunction< GridFunction >( module, "GridFunction" );')
            writer.emit('module.def( "get", [] () { return new GridFunction; } );')
            writer.closePythonModule(name)
            writer.closeNameSpace('FemPy')
            writer.closeNameSpace('Dune')

            writer.close()

            cmake = subprocess.Popen(['cmake', '--build', '../../..', '--target', 'gridfunction'], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, 'gridfunction.so'), os.path.join(compilePath, name + '.so'))

        comm.barrier()
        return importlib.import_module('dune.generated.' + name)
