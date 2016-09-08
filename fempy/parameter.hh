#ifndef DUNE_FEMPY_PARAMETER_HH
#define DUNE_FEMPY_PARAMETER_HH

#include <string>
#include <vector>
#include <utility>

#include <dune/fem/io/parameter.hh>
#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // noParameter
    // -----------

    inline static Fem::ParameterReader noParameter ()
    {
      return Fem::ParameterReader( [] ( const std::string &, const std::string *def ) { return def; } );
    }



    inline static Fem::ParameterReader parameter ( std::vector< std::pair< std::string, std::string > > params )
    {
      typedef std::pair< std::string, std::string > Entry;
      return Fem::ParameterReader( [ params ] ( const std::string &key, const std::string *def ) {
          const auto pos = std::find_if( params.begin(), params.end(), [ &key ] ( const Entry &e ) { return (key == e.first); } );
          return (pos != params.end() ? &pos->second : def);
        } );
    }



    // pyParameter
    // -----------

    inline static Fem::ParameterReader pyParameter ( const pybind11::dict &dict, std::shared_ptr< std::string > tmp )
    {
      return Fem::ParameterReader( [ dict, tmp ] ( const std::string &key, const std::string *def ) {
          pybind11::object value = dict[ key.c_str() ];
          *tmp = (value ? static_cast< std::string >( value.str() ) : Fem::Parameter::getValue( key, *def ));
          return tmp.get();
        } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PARAMETER_HH
