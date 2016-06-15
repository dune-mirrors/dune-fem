#include <config.h>

#include <dune/fem/misc/mpimanager.hh>

#include <dune/fempy/py/base.hh>

#include <dune/fempy/pybind11/pybind11.h>

// VTKDataType
// -----------

PYBIND11_PLUGIN( femmpi )
{
  pybind11::module module( "femmpi" );

  try
  {
    int argc = 0;
    char **argv = nullptr;
    Dune::Fem::MPIManager::initialize( argc, argv );

    // managed to initialize mpi
    typedef Dune::Fem::MPIManager::CollectiveCommunication Comm;

    pybind11::class_< Comm > cc( module, "CollectiveCommunication" );
    cc.def_property_readonly( "rank", &Comm::rank );
    cc.def_property_readonly( "size", &Comm::size );
    cc.def( "barrier", &Comm::barrier );

    module.attr( "comm" ) = pybind11::cast( Dune::Fem::MPIManager::comm() );

    pybind11::enum_< Dune::FemPy::VTKDataType > vtkDataType( module, "DataType" );
    vtkDataType.value( "CellData", Dune::FemPy::VTKDataType::CellData );
    vtkDataType.value( "PointData", Dune::FemPy::VTKDataType::PointData );
    vtkDataType.export_values();
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  return module.ptr();
}
