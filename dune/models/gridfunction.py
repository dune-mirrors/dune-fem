from __future__ import print_function

import importlib
import os
import subprocess
import sys
import types
from dune import comm
from dune.source import SourceWriter

def importFunction(grid, code):
    compilePath = os.path.join(os.path.dirname(__file__), "../generated")

    signature = 'tmp'

    if not isinstance(grid, types.ModuleType):
        grid = grid._module
    name = 'gridfunction_' + signature + '_' + grid._typeHash
    cname = 'GridFunction_' + signature + '_' + grid._typeHash

    if comm.rank == 0:
        if not os.path.isfile(os.path.join(compilePath, name + ".so")):
            writer = SourceWriter(compilePath + '/gridfunction.hh')
            writer.emit(grid._includes)
            writer.emit('')
            writer.emit('#include <dune/corepy/pybind11/pybind11.h>')
            writer.emit('#include <dune/corepy/pybind11/extensions.h>')
            writer.emit('')

            writer.typedef(grid._typeName, 'GridPart')
            writer.typedef("typename GridPart::template Codim< 0 >::EntityType", "EntityType")
           # writer.typedef('Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< double, ' + str(coefficient['dimRange']) + ' > >, GridFunction')
            writer.openStruct(cname, targs=(['class Grid'] + ['class... Coefficients']))

            writer.openConstMethod('bool init', args=['const EntityType &entity'])
            writer.emit('entity_ = &entity;')
            writer.emit('return true;')
            writer.closeConstMethod()

            writer.openConstMethod('const EntityType &entity')
            writer.emit('return *entity_;')
            writer.closeConstMethod()

            writer.openMethod('void evaluate')
            writer.emit(code)
            writer.closeMethod()

            writer.section('private')
            writer.openConstMethod('void initCoefficients', targs=['std::size_t... i'], args=['std::index_sequence< i... >'])
            writer.emit('std::ignore = std::make_tuple( (std::get< i >( coefficients_ ).init( entity() ), i)... );')
            writer.closeConstMethod()
            writer.emit('')
            writer.emit('mutable const EntityType *entity_ = nullptr;')
            writer.emit('mutable std::tuple< Coefficients... > coefficients_;')
            writer.closeStruct()

            writer.emit('')
            writer.typedef(cname + '< GridPart >', 'GridFunction')

            writer.openPythonModule(name)
            writer.emit('')
            writer.emit('// export function class')
            writer.emit('auto cls = pybind11::class_< GridFunction >( module, "GridFunction" );')
            writer.emit('')
            writer.emit('cls.def( "evaluate", [] ( GridFunction &gf ) { return 1; } );')
            writer.emit('module.def( "get", [] () { return new GridFunction(); } );')
            writer.closePythonModule(name)

            writer.close()

            cmake = subprocess.Popen(["cmake", "--build", "../../..", "--target", "gridfunction"], cwd=compilePath)
            cmake.wait()
            os.rename(os.path.join(compilePath, "gridfunction.so"), os.path.join(compilePath, name + ".so"))

        comm.barrier()
        return importlib.import_module("dune.generated." + name)
