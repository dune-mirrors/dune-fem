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

    // noParameter (not used - so can be removed?)
    // -------------------------------------------

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
    namespace detail {

      template <bool logging>
      inline static Fem::ParameterReader pyParameterImpl ( const std::string &rmPrefix,
          std::map<std::string,std::string>&& dict, std::shared_ptr< std::string > tmp,
          const std::string loggingName )
      {
        return Fem::ParameterReader( [ rmPrefix, dict, tmp, loggingName ] ( const std::string &key, const std::string *def ) -> const std::string * {
          // first determine if the `prefix` of the provided key corresponds
          // to the prefix to be removed:
          if (key.compare(0,rmPrefix.size(),rmPrefix) == 0)
          {
            // check the provided map - stripping prefix
            assert(key.size()>rmPrefix.size());
            std::string reducedKey = key.substr(rmPrefix.size(),key.size());
            auto it = dict.find( reducedKey );
            if( it != dict.end() )
            {
              if constexpr (logging)
                Dune::Fem::Parameter::localParameterLog()[loggingName].insert(*it);
              *tmp = it->second;
              return tmp.get();
            }
          }
          else
          {
            auto it = dict.find( key );
            if( it != dict.end() )
            {
              if constexpr (logging)
                Dune::Fem::Parameter::localParameterLog()[loggingName].insert(*it);
              *tmp = it->second;
              return tmp.get();
            }
          }
          // need to check global parameter set
          // the key either does not have the correct prefix or it was not
          // found in the provided map so check the global Parameter container
          if( !Fem::Parameter::exists( key ) )
          {
            if (def == nullptr)
            {
              Dune::Fem::Parameter::localParameterLog()[loggingName].
                       insert({key,""});
              return nullptr; // not found and no default
            }
            else if( *def == Dune::Fem::checkParameterExistsString() )
            {
              return nullptr;
            }
            if constexpr (logging)
              Dune::Fem::Parameter::localParameterLog()[loggingName].
                       insert({key,*def});
            return def;
          }
          if (def)
            Fem::Parameter::get( key, *def, *tmp );
          else
            Fem::Parameter::get( key, *tmp );
          if constexpr (logging)
            Dune::Fem::Parameter::localParameterLog()[loggingName].
                     insert({key,*tmp+" (global)"});
          return tmp.get();
        } );
      }
    } // namespace detail


    inline static Fem::ParameterReader pyParameter ( const std::string &rmPrefix, const pybind11::dict &dict, std::shared_ptr< std::string > tmp )
    {
      typedef std::map< std::string, std::string > DictMap;

      // convert to std::map for enabling find
      DictMap dictMap;
      for( const auto& item : dict )
      {
        dictMap.emplace(
            static_cast< std::string >( pybind11::str(item.first) ),
            static_cast< std::string >( pybind11::str(item.second ) ) );
      }

      if (auto it = dictMap.find( "logging" ); it != dictMap.end() && !it->second.empty() )
        return detail::pyParameterImpl<true>( rmPrefix, std::move(dictMap), tmp, it->second );
      else
        return detail::pyParameterImpl<false>( rmPrefix, std::move(dictMap), tmp, "" );
    }

    inline static Fem::ParameterReader pyParameter ( const pybind11::dict &dict, std::shared_ptr< std::string > tmp )
    {
      return pyParameter("", dict, tmp);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PARAMETER_HH
