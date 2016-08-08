#ifndef DUNE_FEMPY_PARAMETER_HH
#define DUNE_FEMPY_PARAMETER_HH

#include <string>

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

    inline static Fem::ParameterReader pyParameter (const pybind11::dict &dict, std::shared_ptr<std::string> tmp)
    {
      return Fem::ParameterReader( [dict,tmp] ( const std::string &key, const std::string *def )
           { for ( auto entry : dict ) {
               if ( key == static_cast<const std::string&>(entry.first.str()) )
               { *tmp = static_cast<const std::string&>( entry.second.str() ); return tmp.get(); }
             }
             *tmp = Dune::Fem::Parameter::getValue(key,*def); return tmp.get();
           });
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PARAMETER_HH
