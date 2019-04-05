#ifndef DUNE_FEM_IO_PARAMETER_READER_HH
#define DUNE_FEM_IO_PARAMETER_READER_HH

#include <cassert>

#include <functional>
#include <iostream>
#include <string>
#include <utility>

#include <dune/fem/io/parameter/exceptions.hh>
#include <dune/fem/io/parameter/parser.hh>

namespace Dune
{

  namespace Fem
  {

    static inline const std::string& checkParameterExistsString()
    {
      static const std::string defaultKeyForExistCheck("__ParameterReader::check-exists__");
      return defaultKeyForExistCheck;
    }

    // BasicParameterReader
    // --------------------

    template< class Parameter >
    struct BasicParameterReader
    {
      typedef BasicParameterReader<Parameter> ThisType;
      explicit BasicParameterReader ( Parameter parameter = Parameter() )
        : parameter_( std::move( parameter ) )
      {}

      /**
       * \brief check, whether a parameter is defined
       *
       * \param[in]   key    name of the parameter to check
       *
       * \returns \b true, if the parameter is found, \b false otherwise
       */
      bool exists ( const std::string &key ) const
      {
        return static_cast< bool >( parameter_( key, &checkParameterExistsString() ) );
      }

      /**
       * \brief get mandatory parameter
       *
       * \note This method throws an exception, if the parameter cannot be
       *       found.
       *
       * \param[in]   key    name of the parameter to get
       * \param[out]  value  value of the parameter
       */
      template< class T >
      void get ( const std::string &key, T &value ) const
      {
        const std::string *string = parameter_( key, nullptr );
        if( !string )
          DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found." );
        if( !ParameterParser< T >::parse( *string, value ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      /**
       * \brief get optional parameter
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       * \param[in]   key           name of the parameter to get
       * \param[in]   defaultValue  default value for this parameter
       * \param[out]  value         value of the parameter
       */
      template< class T >
      void get ( const std::string &key, const T &defaultValue, T &value ) const
      {
        const std::string defaultString = ParameterParser< T >::toString( defaultValue );
        const std::string *string = parameter_( key, &defaultString );
        assert( string );
        if( !ParameterParser< T >::parse( *string, value ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      /**
       * \brief get optional parameter (special case for string)
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       * \param[in]   key           name of the parameter to get
       * \param[in]   defaultValue  default value for this parameter
       * \param[out]  value         value of the parameter
       */
      void get ( const std::string &key, const char* defaultValue, std::string &value ) const
      {
        const std::string defaultString( defaultValue );
        const std::string *string = parameter_( key, &defaultString );
        assert( string );
        if( !ParameterParser< std::string >::parse( *string, value ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      /**
       * \brief get mandatory parameter
       *
       * \note This method throws an exception, if the parameter cannot be
       *       found.
       *
       * \param[in]   key        name of the parameter to get
       * \param[in]   validator  validator for the parameter value
       * \param[out]  value      value of the parameter
       */
      template< class T, class Validator >
      void getValid ( const std::string &key, const Validator &validator, T &value ) const
      {
        const std::string *string = parameter_( key, nullptr );
        if( !string )
          DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found." );
        if( !ParameterParser< T >::parse( *string, value ) || !validator( value ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      /**
       * \brief get optional parameter
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       * \param[in]   key           name of the parameter to get
       * \param[in]   defaultValue  default value for this parameter
       * \param[in]   validator     validator for the parameter value
       * \param[out]  value         value of the parameter
       */
      template< class T, class Validator >
      void getValid ( const std::string &key, const T &defaultValue, const Validator &validator, T &value ) const
      {
        const std::string defaultString = ParameterParser< T >::toString( defaultValue );
        const std::string *string = parameter_( key, &defaultString );
        assert( string );
        if( !ParameterParser< T >::parse( *string, value ) || !validator( value ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }

      /**
       * \brief get mandatory parameter
       *
       * \note This method throws an exception, if the parameter cannot be
       *       found.
       *
       * \param[in]   key    name of the parameter to get
       *
       * \returns value of the parameter
       */
      template< class T >
      T getValue ( const std::string &key ) const
      {
        T value;
        get( key, value );
        return value;
      }

      /**
       * \brief get optional parameter
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       *  \param[in]   key           name of the parameter to get
       *  \param[in]   defaultValue  default value for this parameter
       *
       *  \returns value of the parameter
       */
      template< class T >
      T getValue ( const std::string &key, const T &defaultValue ) const
      {
        T value = defaultValue;
        get( key, defaultValue, value );
        return value;
      }

      /**
       * \brief get optional parameter
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       *  \param[in]   key           name of the parameter to get
       *  \param[in]   validator     validator for the parameter value
       *
       *  \returns value of the parameter
       */
      template< class T, class Validator >
      T getValidValue ( const std::string &key, const Validator &validator ) const
      {
        T value;
        getValid( key, validator, value );
        return value;
      }

      /**
       * \brief get optional parameter
       *
       * \note This method returns a default value, if the parameter cannot be
       *       found.
       *
       * \param[in]   key           name of the parameter to get
       * \param[in]   defaultValue  default value for this parameter
       * \param[in]   validator     validator for the parameter value
       *
       * \returns value of the parameter
       */
      template< class T, class Validator >
      T getValidValue ( const std::string &key, const T &defaultValue, const Validator &validator ) const
      {
        T value;
        getValid( key, defaultValue, validator, value );
        return value;
      }

      template< int n >
      int getEnum ( const std::string &key, const std::string (&values)[ n ] ) const
      {
        const std::string *string = parameter_( key, nullptr );
        if( !string )
          DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found." );
        return getEnumeration( key, *string, values );
      }

      template< int n >
      int getEnum ( const std::string &key, const std::string (&values)[ n ], int defaultValue ) const
      {
        const std::string *string = parameter_( key, &values[ defaultValue ] );
        return getEnumeration( key, *string, values );
      }

      int getEnum ( const std::string &key, const std::vector<std::string> &values ) const
      {
        const std::string *string = parameter_( key, nullptr );
        if( !string )
          DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found." );
        return getEnumeration( key, *string, values );
      }

      int getEnum ( const std::string &key, const std::vector<std::string> &values, int defaultValue ) const
      {
        const std::string *string = parameter_( key, &values[ defaultValue ] );
        return getEnumeration( key, *string, values );
      }

      ThisType* clone() const { return new ThisType(parameter_); }
      Parameter parameter() { return parameter_; }
      const Parameter parameter() const { return parameter_; }
      void reset() {}

    private:
      template< int n >
      static int getEnumeration ( const std::string &key, const std::string& value, const std::string (&values)[ n ] )
      {
        for( int i = 0; i < n; ++i )
        {
          if( value == values[ i ] )
            return i;
        }

        int j = -1;
        if( !ParameterParser< int >::parse( value, j ) )
          j = -1;
        if( (j < 0) || (j >= n) )
        {
          std::cerr << std::endl << "Parameter '" << key << "' invalid." << std::endl;
          std::cerr << "Valid values are: ";
          for( int i = 0; i < n; ++i )
            std::cerr << values[ i ] << (i < n-1 ? ", " : "");
          std::cerr << std::endl << std::endl;
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
        }
        return j;
      }
      static int getEnumeration ( const std::string &key, const std::string& value, const std::vector<std::string> &values )
      {
        for( unsigned int i = 0; i < values.size(); ++i )
        {
          if( value == values[ i ] )
            return i;
        }

        int j = -1;
        if( !ParameterParser< int >::parse( value, j ) )
          j = -1;
        if( (j < 0) || (j >= (int)values.size()) )
        {
          std::cerr << std::endl << "Parameter '" << key << "' invalid." << std::endl;
          std::cerr << "Valid values are: ";
          for( unsigned int i = 0; i < values.size(); ++i )
            std::cerr << values[ i ] << (i < values.size()-1 ? ", " : "");
          std::cerr << std::endl << std::endl;
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
        }
        return j;
      }

    protected:
      Parameter parameter_;
    };



    // ParameterReader
    // ---------------

    typedef BasicParameterReader< std::function< const std::string *( const std::string &, const std::string * ) > > ParameterReader;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_READER_HH
