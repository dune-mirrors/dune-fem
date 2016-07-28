#include <config.h>

#include <sstream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/corepy/pybind11/pybind11.h>

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
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  {
    pybind11::class_< Dune::Fem::ParameterContainer > param( module, "Parameter" );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::string &filename)
        { Dune::Fem::Parameter::append( filename ); }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::map<std::string,std::string> &dict)
        { for ( auto entery : dict ) Dune::Fem::Parameter::append( entery.first, entery.second ); }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const pybind11::dict &dict)
        { for ( auto entery : dict ) Dune::Fem::Parameter::append( entery.first.str(), entery.second.str() ); }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::string &key, const std::string &val)
        { Dune::Fem::Parameter::append( key, val ); }
    );
    param.def("__str__", [](const Dune::Fem::ParameterContainer &)
        { std::stringstream str; Dune::Fem::Parameter::write( str ); return str.str(); }
    );

    module.attr( "parameter" ) = pybind11::cast( Dune::Fem::Parameter::container() );
  }

  return module.ptr();
}
