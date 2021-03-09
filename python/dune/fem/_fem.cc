#include <config.h>

#include <cmath>
#include <sstream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/threads/threadmanager.hh>

#include <dune/python/pybind11/extensions.h>
#include <dune/python/pybind11/pybind11.h>

PYBIND11_MODULE( _fem, module )
{
  try
  {
    int argc = 0;
    char **argv = nullptr;
    Dune::Fem::MPIManager::initialize( argc, argv );

    int numThreads = 1;
#ifdef USE_SMP_PARALLEL
    {
      const char* nThreads = getenv("OMP_NUM_THREADS");
      if( nThreads )
      {
        numThreads = std::max( int(1), atoi( nThreads ) );
      }
    }
#endif
    Dune::Fem::ThreadManager::setMaxNumberThreads( numThreads );

    if( !pybind11::already_registered< Dune::Fem::MPIManager::CollectiveCommunication >() )
      DUNE_THROW( Dune::Exception, "CollectiveCommunication not registered, yet" );

    module.attr( "comm" ) = pybind11::cast( Dune::Fem::MPIManager::comm() );
  }
  catch ( const std::exception &e )
  {
    std::cout << e.what() << std::endl;
  }

  {
    using pybind11::operator""_a;
    using pybind11::str;

    pybind11::class_< Dune::Fem::ParameterContainer > param( module, "Parameter" );

    param.def( "write", [] ( Dune::Fem::ParameterContainer &self, const std::string &fileName ) {
        std::ofstream file( fileName );
        if( file )
        {
          self.write( file );
          file.close();
        }
      }, "fileName"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &fileName ) {
        self.append( fileName );
      }, "fileName"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::map< std::string, std::string > &entries ) {
        for( auto entry : entries )
          self.append( entry.first, entry.second );
      }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const pybind11::dict &entries ) {
        for ( auto entry : entries )
        {
          std::string s = str( entry.second );
          if (s == "False")
            self.append( str( entry.first ), "false" );
          else if (s == "True")
            self.append( str( entry.first ), "true" );
          else
            self.append( str( entry.first ), s );
        }
      }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
        self.append( key, str( value ) );
      }, "key"_a, "value"_a );

#if 0
    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, int value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );

    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, double value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );
#endif

    param.def( "exists", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
        return self.exists( key );
      }, "key"_a );

    param.def( "__getitem__", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
        if (!self.exists( key ))
          throw pybind11::key_error("key not found in parameter file");
        return self.getValue< std::string >( key );
      } , "key"_a );

    param.def( "__setitem__", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
        self.append( key, str( value ) );
      }, "key"_a, "value"_a );

    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, bool defaultValue ) -> bool {
        if (self.exists( key ))
          return self.getValue< bool >( key );
        else
          self.append( key, defaultValue?"true":"false" );
        return defaultValue;
      }, "key"_a, "defaultValue"_a );
    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, int defaultValue ) {
        if (self.exists( key ))
          return self.getValue< int >( key );
        else
          self.append( key, std::to_string(defaultValue) );
        return defaultValue;
      }, "key"_a, "defaultValue"_a );
    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, double defaultValue ) {
        if (self.exists( key ))
          return self.getValue< double >( key );
        else
          self.append( key, std::to_string(defaultValue) );
        return defaultValue;
      }, "key"_a, "defaultValue"_a );

    param.def( "get", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle defaultValue ) {
        if (self.exists( key ))
          return str( self.getValue< std::string >( key ) );
        else if (!defaultValue.is_none())
          self.append( key, str( defaultValue ) );
        else
          throw pybind11::key_error("key not found in parameter file");
        return str( defaultValue );
      }, "key"_a, pybind11::arg("defaultValue")=pybind11::none() );

    param.def( "__str__", [] ( const Dune::Fem::ParameterContainer &self ) {
        std::stringstream s;
        self.write( s );
        return s.str();
      } );

    module.attr( "parameter" ) = pybind11::cast( Dune::Fem::Parameter::container(),
           pybind11::return_value_policy::reference );
  }
}
