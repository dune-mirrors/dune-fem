#ifndef DUNE_FEMPY_PARAMETER_HH
#define DUNE_FEMPY_PARAMETER_HH

#include <string>
#include <vector>
#include <utility>

#include <dune/fem/io/parameter.hh>
#include <dune/fempy/pybind11/pybind11.hh>

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

    inline static Fem::ParameterReader pyParameter ( const std::string &rmPrefix, const pybind11::dict &dict, std::shared_ptr< std::string > tmp )
    {
      return Fem::ParameterReader( [ rmPrefix, dict, tmp ] ( const std::string &key, const std::string *def ) -> const std::string * {
        // first determine if the `prefix` of the provided key corresponds
        // to the prefix to be removed:
        if (key.compare(0,rmPrefix.size(),rmPrefix) == 0)
        {
          // check the provided map - stripping prefix
          assert(key.size()>rmPrefix.size());
          std::string reducedKey = key.substr(rmPrefix.size(),key.size());
          try {
            pybind11::object value = dict[ reducedKey.c_str() ];
            *tmp = static_cast< std::string >( pybind11::str(value) );
            return tmp.get();
          } catch( ... ) {}
        }
        else
        {
          // check dict without removing prefix
          try {
            pybind11::object value = dict[ key.c_str() ];
            *tmp = static_cast< std::string >( pybind11::str(value) );
            return tmp.get();
          } catch( ... ) {}
        }
        // need to check global parameter set
        // the key either does not have the correct prefix or it was not
        // found in the provided map so check the global Parameter container
        if( !Fem::Parameter::exists( key ) )
        {
          if (def == nullptr)
            return nullptr; // not found and no default
          else if( *def == Dune::Fem::checkParameterExistsString() )
            return nullptr;
          return def;
        }
        if (def)
          Fem::Parameter::get( key, *def, *tmp );
        else
          Fem::Parameter::get( key, *tmp );
        return tmp.get();
      } );
    }
    inline static Fem::ParameterReader pyParameter ( const pybind11::dict &dict, std::shared_ptr< std::string > tmp )
    {
      return pyParameter("", dict, tmp);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PARAMETER_HH
