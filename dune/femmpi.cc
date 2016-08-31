#include <config.h>

#include <sstream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/corepy/pybind11/extensions.h>
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

    if( !pybind11::already_registered< Dune::Fem::MPIManager::CollectiveCommunication >() )
      DUNE_THROW( Dune::Exception, "CollectiveCommunication not registered, yet" );

    module.attr( "comm" ) = pybind11::cast( Dune::Fem::MPIManager::comm() );
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  {
    pybind11::class_< Dune::Fem::ParameterContainer > param( module, "Parameter" );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::string &filename)
        { // int old = Dune::Fem::Parameter::setVerboseRank(0);
          Dune::Fem::Parameter::append( filename );
          // Dune::Fem::Parameter::setVerboseRank( old );
        }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::map<std::string,std::string> &dict)
        { // int old = Dune::Fem::Parameter::setVerboseRank(0);
          for ( auto entery : dict ) Dune::Fem::Parameter::append( entery.first, entery.second );
          // Dune::Fem::Parameter::setVerboseRank( old );
        }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const pybind11::dict &dict)
        { // int old = Dune::Fem::Parameter::setVerboseRank(0);
          for ( auto entery : dict ) Dune::Fem::Parameter::append( entery.first.str(), entery.second.str() );
          // Dune::Fem::Parameter::setVerboseRank( old );
        }
    );
    param.def("append", [](const Dune::Fem::ParameterContainer &, const std::string &key, const pybind11::handle &val)
        { // int old = Dune::Fem::Parameter::setVerboseRank(0);
          Dune::Fem::Parameter::append( key, val.str() );
          // Dune::Fem::Parameter::setVerboseRank( old );
        }
    );
    param.def("__str__", [](const Dune::Fem::ParameterContainer &)
        { std::stringstream str; Dune::Fem::Parameter::write( str ); return str.str(); }
    );

    module.attr( "parameter" ) = pybind11::cast( Dune::Fem::Parameter::container() );
  }

  return module.ptr();
}
