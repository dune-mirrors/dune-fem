#ifndef DUNE_FEM_IO_PARAMETER_CONTAINER_HH
#define DUNE_FEM_IO_PARAMETER_CONTAINER_HH

#include <cassert>
#include <cstddef>

#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <utility>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter/exceptions.hh>
#include <dune/fem/io/parameter/reader.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{

  namespace Fem
  {

    // ParameterContainerData
    // ----------------------

    struct ParameterContainerData
    {
      //! print iteration count and residual information
      static const int solverStatistics   = 1; // this is the new default level
      //! some solver packages offer extended output, use this level for this
      static const int extendedStatistics = 2;
      //! print which parameters have been read, i.e. fem.dofmanager.memoryfactor
      static const int parameterOutput    = 3;
      //! more diagnostics, i.e. about timing and other things
      static const int diagnosticsOutput  = 4;
      //! print debug output at this level
      static const int debugOutput        = 5;

      // default verbosity level used when verbose()
      // without specifying a level is called (old behavior)
      static const int defaultVerbosityLevel = parameterOutput;

      struct Value
      {
        enum ShadowStatus { unresolved, resolved, resolving };

        Value () = default;

        Value ( std::string v, std::string fn ) : value( std::move( v ) ), fileName( std::move( fn ) ) {}

        std::string value, fileName, defaultValue;
        bool used = false, hasDefault = false;
        ShadowStatus shadowStatus = unresolved;
      };

      const std::string *operator() ( const std::string &key, const std::string *defaultValue ) const;

      static std::string trim ( const std::string &s )
      {
        const std::size_t first = s.find_first_not_of( " \t\n" );
        return (first != s.npos ? s.substr( first, s.find_last_not_of( " \t\n" ) + 1 - first ) : std::string());
      }

      std::string resolveEscape ( const std::string &key, std::string &value ) const;
      void resolveShadows ( const std::string &key, Value &val ) const;
      std::string getShadowKey ( const std::string key, const char delimter, std::string &value ) const;

      bool verbose ( const int level = defaultVerbosityLevel ) const
      {
        // return true if verboserank is the current rank and if
        // the activated verbosity level is higher or equal to the given level
        return (verboseRank == MPIManager::rank() && level <= verbosityLevel);
      }

      mutable std::map< std::string, Value > map;
      std::set< std::string > deprecated;
      int verboseRank = 0; // default is to output on rank 0
      int verbosityLevel = 1; // default verbosity level is 1
      bool verbosityLevelPresent = false; // this is parameter was provided
      bool verbosityChangedByVerboseRank = false; // this is true if verboserank was provided, but not verbositylevel
    };



    // ParameterContainer
    // ------------------

    class ParameterContainer
      : public BasicParameterReader< ParameterContainerData >
    {
      typedef ParameterContainerData::Value Value;

      struct DGFBlock;

      static std::string stripComment ( const std::string &line );

      const std::string &insert ( const std::string &key, const std::string &value, bool force );
      bool insert ( const std::string &s, std::queue< std::string > &includes );

      void processFile ( const std::string &filename );
      void processIncludes( std::queue< std::string > &includes );

    public:

      /** \brief cast into ParameterReader */
      operator ParameterReader () const { return ParameterReader( std::ref( parameter_ ) ); }

      /**
       * \brief add parameters from the command line
       *
       * This mehtod adds all parameters (strings containing a colon) in the
       * command line to the container. The parameters are then removed from the
       * command line.
       *
       * \param[in]  argc  number of arguments (as given to main)
       * \param[in]  argv  vector of arguments (as given to main)
       */
      void append ( int &argc, char **argv );

      /**
       * \brief add parameters from a file
       *
       * \param[in]  filename  name of the file containing the parameters
       */
      void append ( const std::string &filename )
      {
        processFile( filename );
      }

      /**
       * \brief add a single parameter to the container
       *
       * \param[in]  key    key of the parameter to add
       * \param[in]  value  value of the parameter to add
       * \param[in]  force  replace parameter, if it exists
       */
      void append ( const std::string &key, const std::string &value, bool force = false )
      {
        if( key != "paramfile" )
        {
          curFileName_ = "program code";
          insert( key, value, force );
        }
        else
          append( value );
      }


      /**
       * \brief A helper function to convert numbers to scientific strings
       *
       * \param[in] value     the value to be converted (needs a << operator)
       */
      template <class T>
      std::string toString( const T& value )
      {
        std::stringstream str;
        str << std::scientific;
        str << value;
        return str.str();
      }

      /**
       * \brief add a single Floating number parameter to the container
       *
       * \param[in]  key    key of the parameter to add
       * \param[in]  value  value of the parameter to add
       * \param[in]  force  replace parameter, if it exists
       */
      template<class NumberType, std::enable_if_t< std::is_floating_point< NumberType >::value || std::is_integral< NumberType >::value, int> = 0 >
      void append ( const std::string &key, NumberType value, bool force = false )
      {
        assert( key != "paramfile" );
        curFileName_ = "program code";
        std::string valueString = toString( value );
        insert( key, valueString, force );
      }

      /**
       * \brief add parameters from a DGF file
       *
       * Parameters can also be read from a DGF file containing a 'FemParameter'
       * block.
       *
       * \param[in]  filename  name of the DGF file containing the parameters
       */
      void appendDGF ( const std::string &filename );

      /** \brief clear all parameters */
      void clear () { parameter_.map.clear(); }

      /** \brief obtain the cached value for fem.verbose */
      bool verbose ( const int level = ParameterContainerData::defaultVerbosityLevel ) const
      {
        return parameter_.verbose( level );
      }

      std::string commonInputPath () const
      {
        return getValue( "fem.prefix.input", std::string( "." ) );
      }

      std::string commonOutputPath () const
      {
        return getValue( "fem.prefix", std::string( "." ) );
      }

      /**
       * \brief write the parameter database to a stream
       *
       * This method writes paramters to the given stream.
       * If the second parameter is true all parameters are written;
       * otherwise only used parameters which do not coincide with the default value
       * are written.
       *
       * \note This method is \b not safe for parallel jobs. Parameters are
       *       written on all ranks.
       *
       * \param[in]  out stream for the parameters.
       * \param[in]  writeAll default is true
       */
      void write ( std::ostream &out, bool writeAll = true ) const;

    private:
      std::string curFileName_;
      int curLineNumber_;
    };



    // ParameterContainer::DGFBlock
    // ----------------------------

    struct ParameterContainer::DGFBlock
      : dgf::BasicBlock
    {
      explicit DGFBlock ( std::istream &in ) : BasicBlock( in, "FemParameter" ) {}

      bool advance () { return getnextline(); }
      std::string getLine () const { return line.str(); }
    };



    // Implementation of ParameterContainerData
    // ----------------------------------------

    inline const std::string *ParameterContainerData::operator() ( const std::string &key, const std::string *defaultValue ) const
    {
      if( deprecated.find( key ) != deprecated.end() )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' deprecated" );

      std::map< std::string, Value >::iterator pos;
      if( defaultValue )
      {
        const std::string& defaultValueStr = *defaultValue;
        // only check existence, do not check default values and the like
        // when the default string has the value of checkParameterExistsString
        // this is to avoid problems with default and non-default parameters
        if( defaultValueStr == checkParameterExistsString() )
        {
          pos = map.find( key );
          if( pos == map.end() )
            return nullptr;
          else
          {
            Value &val = pos->second;
            return &val.value ;
          }
        }

        auto info = map.insert( std::make_pair( key, Value( *defaultValue, "default" ) ) );
        if( info.second && verbose() )
          std::cout << "Adding default: " << key << ": " << *defaultValue << std::endl;
        pos = info.first;
      }
      else
        pos = map.find( key );

      if( pos == map.end() )
        return nullptr;
      Value &val = pos->second;

      if( val.used )
      {
        if( val.hasDefault != static_cast< bool >( defaultValue ) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' used with and without default" );
        if( defaultValue && (val.defaultValue != *defaultValue) )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' used with different default values" );
      }
      else
      {
        val.used = true;
        val.hasDefault = static_cast< bool >( defaultValue );
        if( defaultValue )
          val.defaultValue = *defaultValue;
      }

      resolveShadows( key, val );
      return &val.value;
    }


    inline std::string ParameterContainerData::resolveEscape ( const std::string &key, std::string &value ) const
    {
      if( value.empty() )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' contains trailing '$'." );

      const char escapedChar = value[ 0 ];
      value.replace( 0, 1, "" );

      switch( escapedChar )
      {
      case '$':
      case '%':
      case '#':
        return std::string( "" ) + escapedChar;

      case '(':
        {
          auto pos = map.find( getShadowKey( key, ')', value ) );
          if( pos == map.end() )
            DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found" );
          resolveShadows( pos->first, pos->second );
          return pos->second.value;
        }

      case '[':
        return trim( executeCommand( getShadowKey( key, ']', value ) ) );

      default:
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
      }
    }


    inline void ParameterContainerData::resolveShadows ( const std::string &key, Value &val ) const
    {
      std::string &realValue = val.value;
      if( val.shadowStatus == Value::resolved )
        return;

      if ( val.shadowStatus == Value::resolving )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid, contains infinite loop" );

      val.shadowStatus = Value::resolving;
      std::string realValueHelper;
      realValue.swap( realValueHelper );

      while( !realValueHelper.empty() )
      {
        std::size_t startPoint = realValueHelper.find_first_of( '$' );
        realValue += realValueHelper.substr( 0, startPoint );

        if( startPoint == std::string::npos )
          break;

        realValueHelper.replace( 0, startPoint+1, "" );

        realValue += resolveEscape( key, realValueHelper );
      }
      val.shadowStatus = Value::resolved;
    }


    inline std::string ParameterContainerData::getShadowKey ( const std::string key, const char delimiter, std::string &value ) const
    {
      std::string shadowKey;

      while( true )
      {
        std::size_t startPoint = value.find_first_of( std::string( "$" ) + delimiter );

        if( startPoint == std::string::npos )
          DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );

        shadowKey += value.substr( 0, startPoint );
        const char startChar = value[ startPoint ];

        value.replace( 0, startPoint+1, "" );

        if( startChar == delimiter )
          return shadowKey;
        assert( startChar == '$' );

        shadowKey += resolveEscape( key, value );
      }
    }



    // Implementation of ParameterContainer
    // ------------------------------------

    inline const std::string &ParameterContainer::insert ( const std::string &key, const std::string &value, bool force  = false)
    {
      auto pos = parameter_.map.find( key );
      bool paramExists = ( pos != parameter_.map.end() );
      std::string paramValue;
      if( force && paramExists )
      {
        paramValue = pos->second.value;
        if( paramValue == value )
          return value;
        parameter_.map.erase( key );
      }
      auto info  = parameter_.map.insert( std::make_pair( key, Value( value, curFileName_ ) ) );
      Value &val = info.first->second;
      if( key == "fem.verboserank" )
      {
        ParameterParser< int >::parse( val.value, parameter_.verboseRank );
        if( (parameter_.verboseRank < -1) || (parameter_.verboseRank >= MPIManager::size() ) )
          std::cout << "Warning: Parameter 'fem.verboserank' is neither a " << "valid rank nor -1." << std::endl;

        // Restore default behavior:
        // If fem.verboserank is provided, then we set the verbosityLevel to 3
        // to restore the old behavior. Otherwise the level is 1.
        if( ! parameter_.verbosityLevelPresent &&
            parameter_.verbosityLevel < ParameterContainerData::defaultVerbosityLevel )
        {
          parameter_.verbosityLevel = ParameterContainerData::defaultVerbosityLevel;
          parameter_.verbosityChangedByVerboseRank = true;
        }
      }

      if( key == "fem.verbositylevel" )
      {
        // if verbositylevel is provided undo the changed by verboserank
        if( parameter_.verbosityChangedByVerboseRank )
          parameter_.verbosityLevel = 1;

        parameter_.verbosityLevelPresent = true;

        ParameterParser< int >::parse( val.value, parameter_.verbosityLevel );
        if( (parameter_.verbosityLevel < 0) || (parameter_.verbosityLevel >= 10 ) )
          std::cout << "Warning: Parameter 'fem.verbositylevel' is neither a " << "valid level nor 0." << std::endl;
      }

      if( verbose() )
      {
        std::cout << curFileName_ << "[" << curLineNumber_ << "]: ";
        if( !paramExists )
          std::cout << "Adding " << key << " = " << value << std::endl;
        else if ( !force )
          std::cout << "Ignored " << key << " = " << value << ", using " << val.value << std::endl;
        else
          std::cout << "Replacing " << key << " = " << paramValue << " by " << value << std::endl;
      }

      return force ? value : val.value;
    }


    inline std::string ParameterContainer::stripComment ( const std::string &line )
    {
      std::size_t size = line.size();
      std::size_t end = line.find_first_of ( "%#$" );

      while( (end != std::string::npos) && (line[end] =='$') )
      {
        if( end+2 < size )
          end = line.find_first_of ( "%#$", end+2 );
        else
          end = std::string::npos;
      }

      return ParameterContainerData::trim( line.substr( 0, end ) );
    }


    inline bool ParameterContainer::insert ( const std::string &s, std::queue< std::string > &includes )
    {
      const std::size_t size = s.size();

      std::size_t key_start = 0;
      for( ; key_start < size; ++key_start )
      {
        if( (s[ key_start ] != ' ') && (s[ key_start ] != '\t') )
          break;
      }

      std::size_t key_end = key_start;
      for( ; key_end < size; ++key_end )
      {
        const char &c = s[ key_end ];
        if( (c == ' ') || (c == '\t') || (c == ':') )
          break;
      }

      std::size_t value_start = key_end;
      for( ; value_start < size ; ++value_start )
      {
        if( s[ value_start ] == ':' )
          break;
      }
      ++value_start;

      for( ; value_start < size; ++value_start )
      {
        if( (s[ value_start ] != ' ') && (s[ value_start ] != '\t') )
          break;
      }

      std::size_t value_end = value_start;
      for( std::size_t i = 0; i < size; ++i )
      {
        if( (s[ i ] != ' ') && (s[ i ] != '\t') )
          value_end = i+1;
      }

      if( value_start >= size )
        return false;

      std::string key = s.substr( key_start, key_end - key_start );
      std::string value = s.substr( value_start, value_end - value_start );

      if( key == "paramfile" )
        includes.push( commonInputPath() + "/" + value );
      else if( key == "deprecated" )
        parameter_.deprecated.insert( value );
      else
        insert( key, value );
      return true;
    }


    inline void ParameterContainer::processFile ( const std::string &filename )
    {
      if( verbose() )
        std::cout << "Parameter: Processing '" << filename << "'..." << std::endl;

      std::ifstream file( filename );
      if( !file.is_open() )
      {
        std::cerr << "Warning: Unable to read parameter file '" << filename << "'" << std::endl;
        return;
      }

      curFileName_ = filename;
      curLineNumber_ = 0;
      std::queue< std::string > includes;

      while( !file.eof() )
      {
        std::string line;
        std::getline( file, line );
        curLineNumber_++;
        line = stripComment( line );
        if( !line.empty() )
          insert( line, includes );
      }
      file.close();

      processIncludes( includes );
    }


    inline void ParameterContainer::processIncludes( std::queue< std::string > &includes )
    {
      while( !includes.empty() )
      {
        Value val;
        val.value = includes.front();
        includes.pop();
        parameter_.resolveShadows( "paramfile", val );
        processFile( val.value );
      }
    }


    inline void ParameterContainer::append ( int &argc, char **argv )
    {
      std::queue< std::string > includes;
      curFileName_ = "program arguments";
      curLineNumber_ = 0;
      for( int i = 1 ; i < argc; ++i )
      {
        ++curLineNumber_;
        if( !insert( std::string( argv[ i ] ), includes ) )
          continue;

        std::copy( argv + (i+1), argv + argc, argv + i );
        --i;
        --argc;
      }

      processIncludes( includes );
    }


    inline void ParameterContainer::appendDGF ( const std::string &filename )
    {
      if( verbose() )
        std::cout << "Parameter: Processing DGF '" << filename << "'..." << std::endl;

      std::ifstream file( filename );
      if( !file.is_open() )
      {
        std::cerr << "Warning: Unable to read DGF file '" << filename << "'" << std::endl;
        return;
      }

      if( !DuneGridFormatParser::isDuneGridFormat( file ) )
        return;

      DGFBlock block( file );
      if( !block.isactive() )
        return;

      curFileName_ = filename;
      curLineNumber_ = 0;
      std::queue< std::string > includes;

      while( block.advance() )
      {
        ++curLineNumber_;
        const std::string line = stripComment( block.getLine() );
        if( !line.empty() )
          insert( line, includes );
      }

      processIncludes( includes );
    }


    inline void ParameterContainer::write ( std::ostream &out, bool writeAll ) const
    {
      std::map< std::string, std::map<std::string, std::string> > writeMap;
      for( const auto &param : parameter_.map )
      {
        const Value &val = param.second;
        if( writeAll || !val.hasDefault || (val.used && (val.value != val.defaultValue)) )
          writeMap[ val.fileName ][ (val.used ? "": "# " ) + param.first ] = val.value;
      }

      for( const auto &source : writeMap )
      {
        out << "# from " << source.first << std::endl;
        for( const auto &param : source.second )
          out << param.first << ": " << param.second << std::endl;
        out << std::endl;
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_PARAMETER_CONTAINER_HH
