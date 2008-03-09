#ifndef DUNE_FEM_PARAMETER_HH
#define DUNE_FEM_PARAMETER_HH

#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/exceptions.hh>

namespace Dune
{

  class ParameterNotFound
  : public Exception
  {};

  class ParameterInvalid
  : public Exception
  {};



  template< class T >
  class ValidateNotGreater
  {
  protected:
    const T threshold_;

  public:
    inline ValidateNotGreater ( const T &threshold )
    : threshold_( threshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value <= threshold_;
    }
  };
  template< class T >
  class ValidateNotLess
  {
  protected:
    const T threshold_;

  public:
    inline ValidateNotLess ( const T &threshold )
    : threshold_( threshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value >= threshold_;
    }
  };
  template< class T >
  class ValidateIsLess
  {
  protected:
    const T threshold_;

  public:
    inline ValidateIsLess ( const T &threshold )
    : threshold_( threshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value < threshold_;
    }
  };
  template< class T >
  class ValidateIsGreater
  {
  protected:
    const T threshold_;

  public:
    inline ValidateIsGreater ( const T &threshold )
    : threshold_( threshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return value > threshold_;
    }
  };
  template< class T >
  class ValidateInOpenInterval
  {
  protected:
    const T lThreshold_,rThreshold_;

  public:
    inline ValidateInOpenInterval ( const T &lThreshold, const T &rThreshold )
    : lThreshold_( lThreshold ),
      rThreshold_( rThreshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return (value > lThreshold_ && value < rThreshold_); 
    }
  };
  template< class T >
  class ValidateInClosedInterval
  {
  protected:
    const T lThreshold_,rThreshold_;

  public:
    inline ValidateInClosedInterval ( const T &lThreshold, const T &rThreshold )
    : lThreshold_( lThreshold ),
      rThreshold_( rThreshold )
    {}

    inline bool operator() ( const T &value ) const
    {
      return (value >= lThreshold_ && value <= rThreshold_); 
    }
  };



  /** \class Parameter
   *  \brief Container for User Specified Parameters
   *
   *  The class Parameter provides parameters collected from given parameter
   *  files or the command line in a unified and easy to use way.
   *
   *  This class adheres to the singleton concept, i.e., all methods are static
   *  and internally use a single instance of this object to store all data.
   */
  class Parameter
  {
    typedef Parameter ThisType;

  private:
    std :: map< std :: string, std :: string > params_;
    bool verbose_;
    
  private:
    inline Parameter ()
    : verbose_( false )
    {}

    inline Parameter ( const ThisType & );
    inline ThisType &operator= ( const ThisType & );

    inline static ThisType &instance ()
    {
      static ThisType theInstance;
      return theInstance;
    }

    inline const std :: string &map ( const std :: string &key )
    {
      if( params_.find( key ) == params_.end() )
      {
        std :: ostringstream message;
        message << "Parameter '" << key << "' not found.";
        DUNE_THROW( ParameterNotFound, message.str() );
      }
      return params_[ key ];
    }

    inline const std :: string &map ( const std :: string &key,
                                      const std :: string &value )
    {
      if( params_.find( key ) == params_.end() )
      {
        if( verbose_ )
          std :: cout << "Parameter: Adding " << key << " = " << value
                      << std :: endl;
        params_[ key ] = value;
      }
      return params_[ key ];
    }

    inline void insert ( const std :: string s,
                         std :: queue< std :: string > &includes )
    {
      const unsigned int size = s.size();

      unsigned int key_start = 0;
      for( ; key_start < size; ++key_start )
      {
        if( (s[ key_start ] != ' ') && (s[ key_start ] != '\t') )
          break;
      }

      unsigned int key_end = key_start;
      for( ; key_end < size; ++key_end )
      {
        const char &c = s[ key_end ];
        if( (c == ' ') || (c == '\t') || (c == ':') )
          break;
      }

      unsigned int value_start = key_end;
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

      if( value_start >= size )
        return;

      std :: string key = s.substr( key_start, key_end - key_start );
      std :: string value = s.substr( value_start, size - value_start );

      if( key != "paramfile" )
      {
        const std :: string &actual_value = map( key, value );
        if( key == "fem.verbose" )
          parse( actual_value, verbose_ );
      }
      else
        includes.push( value );
    }

    template< class T >
    inline static void parse ( const std :: string &s, T &value )
    {
      std :: istringstream in( s );
      in >> value;
    }

    inline void processFile ( const std :: string &filename )
    {
      if( verbose_ )
        std :: cout << "Parameter: Processing '" << filename << "'..."
                     << std :: endl;
      
      std :: ifstream file( filename.c_str() );
      if( file.is_open() )
      {
        std :: queue< std :: string > includes;
        while( !file.eof() )
        {
          std :: string line;
          std :: getline( file, line );
          if( line.size() == 0 )
            continue;

          if( (line[ 0 ] != '%') && (line[ 0 ] != '#') )
            insert( line, includes );
        }
        file.close();

        processIncludes( includes );
      }
      else
        std :: cerr << "Warning: Unable to read parameter file '"
                    << filename << "'" << std :: endl;
    }

    inline void processIncludes( std :: queue< std :: string > &includes )
    {
      while( !includes.empty() )
      {
        std :: string filename = includes.front();
        includes.pop();
        processFile( filename );
      }
    }

    inline void replace ( const std :: string &key,
                          const std :: string &value )
    {
      if( verbose_ )
        std :: cout << "Parameter: Replacing " << key << " = " << value
                    << std :: endl;
      params_[ key ] = value;
    }

  public:
    /** \brief add parameters from a file to the container
     * 
     * \param[in]  filename  name of the file containing the parameters
     */
    inline static void append ( const std :: string &filename )
    {
      instance().processFile( filename );
    }

    /** \brief add parameters from the command line to the container
     * 
     * \param[in]  argc  number of arguments (as given to main)
     * \param[in]  argv  vector of arguments (as given to main)
     */
    inline static void append ( int argc, char **argv )
    {
      std :: queue< std :: string > includes;
      for( int i = 1; i < argc; ++i )
        instance().insert( std :: string( argv[ i ] ), includes );
      instance().processIncludes( includes );
    }

    /** \brief get a mandatory parameter from the container
     *
     *  \note This method throws an exception, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key    name of the parameter to get
     *  \param[out]  value  value of the parameter
     */
    template< class T >
    inline static void get ( const std :: string &key,
                             T &value )
    {
      parse( instance().map( key ), value );
    }
    
    /** \brief get an optional parameter from the container
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   defaultValue  default value for this parameter
     *  \param[out]  value         value of the parameter
     */
    template< class T >
    inline static void get ( const std :: string &key,
                             const T &defaultValue,
                             T &value )
    {
      std :: ostringstream out;
      out << defaultValue;
      parse( instance().map( key, out.str() ), value );
    }

    /** \brief get a mandatory parameter from the container
     *
     *  \note This method throws an exception, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key        name of the parameter to get
     *  \param[in]   validator  validator for the parameter value
     *  \param[out]  value      value of the parameter
     */
    template< class T, class Validator >
    inline static void getValid ( const std :: string &key,
                                  const Validator &validator,
                                  T &value )
    {
      parse( instance().map( key ), value );
      if( !validator( value ) )
      {
        std :: ostringstream message;
        message << "Parameter '" << key << "' invalid.";
        DUNE_THROW( ParameterInvalid, message.str() );
      }
    }

    /** \brief get an optional parameter from the container
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   defaultValue  default value for this parameter
     *  \param[in]   validator     validator for the parameter value
     *  \param[out]  value         value of the parameter
     */
    template< class T, class Validator >
    inline static void getValid ( const std :: string &key,
                                  const T &defaultValue,
                                  const Validator &validator,
                                  T &value )
    {
      std :: ostringstream out;
      out << defaultValue;
      parse( instance().map( key, out.str() ), value );
      if( !validator( value ) )
      {
        std :: cerr << "Warning: Parameter '" << key << "' is invalid."
                    << std :: endl;
        instance().replace( key, out.str() );
        parse( out.str(), value );
      }
    }
    
    /** \brief get a mandatory parameter from the container
     *
     *  \note This method throws an exception, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key    name of the parameter to get
     *  
     *  \returns value of the parameter
     */
    template< class T >
    inline static T getValue ( const std :: string &key )
    {
      T value;
      get( key, value );
      return value;
    }
    
    /** \brief get an optional parameter from the container
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   defaultValue  default value for this parameter
     *
     *  \returns value of the parameter
     */
    template< class T >
    inline static T getValue ( const std :: string &key,
                               const T &defaultValue )
    {
      T value;
      get( key, defaultValue, value );
      return value;
    }
    
    /** \brief get an optional parameter from the container
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   validator     validator for the parameter value
     *
     *  \returns value of the parameter
     */
    template< class T, class Validator >
    inline static T getValidValue ( const std :: string &key,
                                    const Validator &validator )
    {
      T value;
      getValid( key, validator, value );
      return value;
    }


    /** \brief get an optional parameter from the container
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   defaultValue  default value for this parameter
     *  \param[in]   validator     validator for the parameter value
     *
     *  \returns value of the parameter
     */
    template< class T, class Validator >
    inline static T getValidValue ( const std :: string &key,
                                    const T &defaultValue,
                                    const Validator &validator )
    {
      T value;
      getValid( key, defaultValue, validator, value );
      return value;
    }

    /** \brief obtain the value for fem.prefix
     */
    inline static const std :: string &prefix ()
    {
      return instance().map( "fem.prefix", "." );
    }

    /** \brief obtain the cached value for fem.verbose
     */
    inline static bool verbose ()
    {
      return instance().verbose_;
    }
  
    /** \brief write all parameters into a file
     *
     *  \param[in]  filename  name of the file to store the parameters in
     */
    inline static void write ( const std :: string &filename )
    {
      typedef std :: map< std :: string, std :: string > :: iterator iterator;

      std :: string fullname( prefix() );
      fullname += "/";
      fullname += filename;

      std :: ofstream file( fullname.c_str() );
      if( !file.is_open() )
      {
        std :: cerr << "Warning: Unable to write parameter file '"
                    << filename << "'" << std :: endl;
      }

      const iterator end = instance().params_.end();
      for( iterator it = instance().params_.begin(); it != end; ++it )
        file << it->first << ": " << it->second << std :: endl;

      file.close();
    }
  };
  
}

#endif
