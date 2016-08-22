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

    signature = 'tmp'

    dimRange = grid.dimensionworld
    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'gridfunction_' + signature + '_' + grid._typeHash
    cname = 'GridFunction_' + signature + '_' + grid._typeHash

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

            writer.typedef(grid._typeName, 'GridPart')
            writer.emit('static const int dimRange = ' + str(dimRange) + ';')
            writer.emit('static const int dimDomain = GridPart::dimensionworld;')
            writer.typedef('typename Dune::Fem::FunctionSpace< double, double, dimDomain, dimRange >', 'FunctionSpaceType')
            writer.typedef('typename GridPart::template Codim< 0 >::EntityType', 'EntityType')
            writer.typedef('typename FunctionSpaceType::DomainType', 'DomainType')
            writer.typedef('typename EntityType::Geometry::LocalCoordinate', 'LocalCoordinateType')
            writer.typedef('typename FunctionSpaceType::RangeType', 'RangeType')
            writer.typedef('typename FunctionSpaceType::JacobianRangeType', 'JacobianRangeType')
            writer.typedef('typename FunctionSpaceType::HessianRangeType', 'HessianRangeType')
            writer.openStruct(cname, targs=(['class GridPart'] + ['class RangeType'] + ['class... Coefficients']))

            writer.openConstMethod('bool init', args=['const EntityType &entity'])
            writer.emit('entity_ = &entity;')
            writer.emit('return true;')
            writer.closeConstMethod()

            writer.openConstMethod('const EntityType &entity')
            writer.emit('return *entity_;')
            writer.closeConstMethod()

            writer.openConstMethod('void evaluate', args=['const DomainType &x', 'RangeType &value'])
            writer.emit(code)
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
            writer.typedef(cname + '< GridPart, RangeType >', 'GridFunction')

            writer.openNameSpace('Dune')
            writer.openNameSpace('FemPy')
            writer.openPythonModule(name)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('')
            #writer.emit('auto cls = pybind11::class_< GridFunction >( module, "GridFunction");')
            writer.emit('pybind11::class_< GridFunction > registerGridFunction( module );')
            writer.closePythonModule(name)
            writer.closeNameSpace('FemPy')
            writer.closeNameSpace('Dune')

            #writer.emit('cls.def_property_readonly( "dimRange", [] ( GridFunction & ) -> int { return RangeType::dimension; } );')
            #writer.emit('cls.def( "evaluate", [] ( const GridFunction &gf, const LocalCoordinateType &x ) {')
            #writer.emit('RangeType value;')
            #writer.emit('gf.evaluate( x, value );')
            #writer.emit('return value; } );')
            #writer.emit('module.def( "get", [] () { return new GridFunction(){}; } );')

            writer.close()

            cmake = subprocess.Popen(['cmake', '--build', '../../..', '--target', 'gridfunction'], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, 'gridfunction.so'), os.path.join(compilePath, name + '.so'))

        comm.barrier()
        return importlib.import_module('dune.generated.' + name)
