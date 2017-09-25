#include <config.h>

#include <sstream>

#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/io/parameter.hh>

#include <dune/corepy/pybind11/extensions.h>
#include <dune/corepy/pybind11/pybind11.h>

PYBIND11_MODULE( _fem, module )
{
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
    using pybind11::operator""_a;
    using pybind11::str;

    pybind11::class_< Dune::Fem::ParameterContainer > param( module, "Parameter" );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &fileName ) {
        self.append( fileName );
      }, "fileName"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::map< std::string, std::string > &entries ) {
        for( auto entry : entries )
          self.append( entry.first, entry.second );
      }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const pybind11::dict &entries ) {
        for ( auto entry : entries )
          self.append( str( entry.first ), str( entry.second ) );
      }, "entries"_a );

    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
        self.append( key, str( value ) );
      }, "key"_a, "value"_a );

    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, int value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );

    // do we really need this one?
    param.def( "append", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, double value ) {
        self.append( key, std::to_string( value ) );
      }, "key"_a, "value"_a );

    param.def( "exists", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
        return self.exists( key );
      }, "key"_a );

    param.def( "__getitem__", [] ( const Dune::Fem::ParameterContainer &self, const std::string &key ) {
        return self.getValue< std::string >( key );
      } , "key"_a );

    param.def( "__setitem__", [] ( Dune::Fem::ParameterContainer &self, const std::string &key, pybind11::handle value ) {
        self.append( key, str( value ) );
      }, "key"_a, "value"_a );

    param.def( "__str__", [] ( const Dune::Fem::ParameterContainer &self ) {
        std::stringstream s;
        self.write( s );
        return s.str();
      } );

    module.attr( "parameter" ) = pybind11::cast( Dune::Fem::Parameter::container(),
           pybind11::return_value_policy::reference );
  }
}
