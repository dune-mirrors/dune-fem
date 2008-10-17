#ifndef DUNE_FEM_PARAMETER_HH
#define DUNE_FEM_PARAMETER_HH

#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/validator.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
  
  /** \addtogroup Parameter
   *
   *  Handling Parameters, i.e., values that can be set after compilation, in
   *  dune-fem is extremely easy. Just add
   *  \code
   *  Dune :: MPIManager :: initialize( argc, argv );
   *  Dune :: Parameter :: append( argc, argv );
   *  \endcode
   *  at the head of your main function. Parameters are strings of the
   *  format "key: value". Any command line argument containing a colon
   *  will thus be interpreted as a parameter.
   *
   *  \note All parameter keys are case sensitive.
   *  \note Found parameters are removed from the command line.
   *
   *  There are 8 static methods in Parameter to obtain the value of a
   *  parameter. They are divided by the following criteria:
   *  - \b Default \b Value: The methods taking a default value will return the
   *    default, if the parameter has not been specified by the user.
   *    Additionally, they will add "\e key : \e value" to the database. The
   *    methods not taking a default value will throw an exception, if the
   *    parameter could not be found in the database.
   *  - \b Return \b Value: For convenience, there is always a method (called
   *    getValue) returning the value of the parameter. If you do not want to
   *    rely on return value optimization, use the method (called get) taking
   *    a reference to the return variable as an argument.
   *  - \b Validation: It is often necessary to make sure the value satisfies
   *    certain constraints. For this purpose, some methods take a validator
   *    as an argument.
   *  .
   *
   *  Of course, you don't have to pass every parameter on the command line.
   *  They can also be gathered in files. Parameter provides a kind of include
   *  mechanism. Whenever a parameter with key "paramfile" is encountered, the
   *  value is interpreted as a paramter file to include.
   *
   *  If a parameter is defined multiply, the first definition is added to the
   *  database. All later definitions are ignored. Therefore it is important to
   *  know the exact behaviour of the "paramfile" parameter:
   *  - All parameters in the current file (or the command line) are added
   *    first.
   *  - If there were includes, they are processed in the order of appearance.
   *  - Should an included file have includes, they are added depth-first,
   *    i.e., The includes of one included file are parsed down to the last
   *    file included, before the includes of the next included file are
   *    considered.
   *  .
   *
   *  All parameter names defined by dune-fem should conform to the following
   *  naming convention:
   *  \code
   *  fem.<group>.<parameter>
   *  \endcode
   *  The group name can be omitted if necessary.
   *  
   *  An example is the parameter 
   *  \code 
   *  fem.verboserank
   *  \endcode
   *  This can beused throughout the the program; by calling:
   *  \code
   *  Parameter :: verbose()
   *  \endcode
   *  If verbose is set, information concerning the parameters read will
   *  be output to stdout.
   *
   *  Here is an example usage:
   *  \code
   *  #include <dune/fem/io/parameter.hh>
   *  
   *  using Dune :: Parameter;
   *
   *  int globalFlag;
   *
   *  int main ( int argc, char **argv )
   *  {
   *    Dune :: MPIManager :: initialize( argc, argv );
   *    Parameter :: append( argc, argv );
   *
   *    // get parameters
   *    double startTime = Parameter :: getValue< double >( "starttime", 0.0 );
   *    Dune :: ValidateGreater< double > validator( startTime );
   *    double endTime = Parameter :: getValidValue< double >( "endtime", validator );
   *    Parameter :: get( "flag", 0, globalFlag );
   *    
   *    if( Parameter :: verbose() )
   *      std :: cout << "Computing from " << startTime << " to " << endTime << std::endl;
   *    
   *    // ...
   *
   *    std :: ofstream results( (Parameter :: outputPrefix() + "/results").c_str() );
   *
   *    // ...
   *
   *    Parameter :: write( "parameter.log" );
   *  }
   *  \endcode
   */

  class ParameterNotFound
  : public Exception
  {};

  class ParameterInvalid
  : public Exception
  {};



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
    typedef std :: map< std :: string, std :: string > ParameterMapType;

  private:
    std :: string curFileName_;
    int curLineNumber_;
    ParameterMapType params_;
    int verboseRank_;
    
  private:
    inline Parameter ()
    : verboseRank_( -1 )
    {}

    // Prohibit copying and assignment
    inline Parameter ( const ThisType & );
    inline ThisType &operator= ( const ThisType & );

    inline static ThisType &instance ()
    {
      static ThisType theInstance;
      return theInstance;
    }

    inline const std :: string *find ( const std :: string &key );

    inline const std :: string &map ( const std :: string &key );

    inline const std :: string &map ( const std :: string &key,
                                      const std :: string &value,
                                      bool verbFound = false );

    inline bool insert ( const std :: string &s,
                         std :: queue< std :: string > &includes );

    template< class T >
    inline static void parse ( const std :: string &s, T &value );

    inline void processFile ( const std :: string &filename );

    inline void processIncludes( std :: queue< std :: string > &includes );

    inline void replace ( const std :: string &key,
                          const std :: string &value );

  public:
    /** \brief add parameters from a file to the container
     * 
     * \param[in]  filename  name of the file containing the parameters
     */
    inline static void append ( const std :: string &filename )
    {
      instance().processFile( filename );
    }

    /** \brief clear all parameters
     */
    inline static void clear() {
      instance().params_.clear();
    }
   
    template <class T>
    inline static void replaceKey ( const std :: string& key, const T& value)
    {
      std::stringstream valueStr;
      valueStr << value;
      instance().replace( key,valueStr.str() );
    }

    /** \brief add parameters from the command line
     *
     *  This mehtod adds all parameters (strings containing a colon) in the
     *  command line to the container. The parameters are then removed from the
     *  command line.
     * 
     * \param[in]  argc  number of arguments (as given to main)
     * \param[in]  argv  vector of arguments (as given to main)
     */
    inline static void append ( int &argc, char **argv );
    
    /** \brief find out, whether a parameter is defined in the container
     *
     *  \param[in]   key    name of the parameter to check
     *
     *  \returns \b true, if the parameter is found in the container,
     *           \b false otherwise
     */
    inline static bool exists ( const std :: string &key )
    {
      return (instance().find( key ) != 0);
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
                             T &value );

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
                                  T &value );

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
                                  T &value );
    
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
    static T getValue ( const std :: string &key )
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
    static T getValue ( const std :: string &key,
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
    static T getValidValue ( const std :: string &key,
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
    static T getValidValue ( const std :: string &key,
                             const T &defaultValue,
                             const Validator &validator )
    {
      T value;
      getValid( key, defaultValue, validator, value );
      return value;
    }

    /** \brief obtain common output path
     *
     *  For parallel jobs you need two different output paths:
     *  - Data common to all processes should be written to commonOutputPath().
     *    Only one process (i.e., rank 0) should write this data.
     *  - Data unique to each process should be written to outputPath().
     *    Each process should write this data.
     *  .
     * 
     *  \returns value of parameter 'fem.prefix', which defaults to '.'.
     */
    static std :: string commonOutputPath ()
    {
      return instance().map( "fem.prefix", "." );
    }

    /** \brief obtain unique output path for this process
     *
     *  For parallel jobs you need two different output paths:
     *  - Data common to all processes should be written to commonOutputPath().
     *    Only one process (i.e., rank 0) should write this data.
     *  - Data unique to each process should be written to outputPath().
     *    Each process should write this data.
     *  .
     *
     *  \returns \<prefix\>/p\<rank\>, where
     *  - \<prefix\> denotes the value of 'fem.prefix', which defaults to '.',
     *  - \<rank\> denots the this processes rank.
     *  .
     */
    inline static std :: string outputPath ();

    /** \brief obtain the value for fem.prefix defaults to '.'
     */
    DUNE_DEPRECATED static const std :: string &prefix ()
    {
      return instance().map( "fem.prefix", "." );
    }

    /** \brief obtain the cached value for fem.verbose
     */
    static bool verbose ()
    {
      return (instance().verboseRank_ == MPIManager :: rank());
    }
  
    /** \brief write the parameter database to a file
     *
     *  This method writes all paramters in the database to the given file.
     *  The parameters are stored in alphabetical order. Includes are not used.
     *
     *  \note This method is safe for parallel jobs. Parameters are only written
     *        on rank 0.
     *
     *  \param[in]  filename  name of the file to store the parameters in; prefix() is used.
     */
    inline static void write ( const std :: string &filename );

    /** \brief write the parameter database to a stream
     *
     *  This method writes all paramters in the database to the given stream.
     *
     *  \note This method is \b not safe for parallel jobs. Parameters are
     *        written on all ranks.
     *
     *  \param[in]  out stream for the parameters.
     */
    inline static void write ( std :: ostream &out );
  };



  // Private Methods
  // ---------------

  inline const std :: string *
  Parameter :: find ( const std :: string &key )
  {
    ParameterMapType :: iterator it = params_.find( key );
    return (it != params_.end()) ? &(it->second) : 0;
  }

  inline const std :: string &
  Parameter :: map ( const std :: string &key )
  {
    const std :: string *value = find( key );
    if( value != 0 )
      return *value;

    std :: ostringstream message;
    message << "Parameter '" << key << "' not found.";
    DUNE_THROW( ParameterNotFound, message.str() );
  }

  inline const std :: string &
  Parameter :: map ( const std :: string &key,
                     const std :: string &value,
                     bool verbFound )
  {
    std :: pair< ParameterMapType :: iterator, bool > info
      = params_.insert( std :: make_pair( key, value ) );

    if( info.second )
    {
      if( verbose() )
      {
        std :: cout << curFileName_ << "[" << curLineNumber_ << "]: ";
        std :: cout << "Adding " << key << " = " << value << std :: endl;
      }
    }
    else
    {
      if( verbFound && verbose() )
      {
        std :: cout << curFileName_ << "[" << curLineNumber_ << "]: ";
        std :: cout << "Ignored " << key << " = " << value
                    << ", using " << info.first->second << std::endl;
      }
    }
    return info.first->second;
  }

  inline bool
  Parameter :: insert ( const std :: string &s,
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
      return false;

    std :: string key = s.substr( key_start, key_end - key_start );
    std :: string value = s.substr( value_start, size - value_start );

    if( key != "paramfile" )
    {
      const std :: string &actual_value = map( key, value, true );
      if( key == "fem.verbose" )
      {
        bool verbose;
        parse( actual_value, verbose );
        verboseRank_ = (verbose ? MPIManager :: rank() : -1);
        std :: cout << "Warning: Using deprecated parameter 'fem.verbose'; "
                    << "use 'fem.verboserank' instead." << std :: endl;
      }
      if( key == "fem.verboserank" )
      {
        parse( actual_value, verboseRank_ );
        if( (verboseRank_ < -1) || (verboseRank_ >= MPIManager :: size() ) )
        {
          std :: cout << "Warning: Parameter 'fem.verboserank' is neither a "
                      << "valid rank nor -1." << std :: endl;
        }
      }
    }
    else
      includes.push( value );
    return true;
  }

  template< class T >
  inline void
  Parameter :: parse ( const std :: string &s, T &value )
  {
    std :: istringstream in( s );
    in >> value;
  }

  inline void
  Parameter :: processFile ( const std :: string &filename )
  {
    if( verbose() )
      std :: cout << "Parameter: Processing '" << filename << "'..."
                   << std :: endl;
    
    std :: ifstream file( filename.c_str() );
    if( file.is_open() )
    {
      curFileName_ = filename;
      curLineNumber_ = 0;
      std :: queue< std :: string > includes;
      while( !file.eof() )
      {
        std :: string line;
        std :: getline( file, line );
        curLineNumber_++;
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

  inline void
  Parameter :: processIncludes( std :: queue< std :: string > &includes )
  {
    while( !includes.empty() )
    {
      std :: string filename = includes.front();
      includes.pop();
      processFile( filename );
    }
  }

  inline void
  Parameter :: replace ( const std :: string &key,
                         const std :: string &value )
  {
    if( verbose() )
      std :: cout << "Parameter: Replacing " << key << " = " << value
                  << std :: endl;
    params_[ key ] = value;
  }



  // Public methods
  // --------------

  inline void
  Parameter :: append ( int &argc, char **argv )
  {
    std :: queue< std :: string > includes;
    instance().curFileName_ = "programm arguments";
    int &i = instance().curLineNumber_;
    for( i = 1 ; i < argc; ++i )
    {
      if( !instance().insert( std :: string( argv[ i ] ), includes ) )
        continue;

      for( int j = i+1; j < argc; ++j )
        argv[ j-1 ] = argv[ j ];
      --i;
      --argc;
    }
    instance().processIncludes( includes );
  }

  template< class T >
  inline void
  Parameter :: get ( const std :: string &key,
                     const T &defaultValue,
                     T &value )
  {
    instance().curFileName_ = "using default";
    instance().curLineNumber_ = 0;
    std :: ostringstream out;
    out << defaultValue;
    parse( instance().map( key, out.str() ), value );
  }

  template< class T, class Validator >
  inline void
  Parameter :: getValid ( const std :: string &key,
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

  template< class T, class Validator >
  inline void
  Parameter :: getValid ( const std :: string &key,
                          const T &defaultValue,
                          const Validator &validator,
                          T &value )
  {
    instance().curFileName_ = "using default";
    instance().curLineNumber_ = 0;
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
  
  inline std :: string
  Parameter :: outputPath ()
  {
    const std :: string &prefix = instance().map( "fem.prefix", "." );
    std :: ostringstream path;
    if( !prefix.empty() )
    {
      path << prefix;
      if( prefix[ prefix.length()-1 ] != '/' )
        path << '/';
    }
    else
      path << "./";
    path << 'p' << MPIManager :: rank();
    return path.str();
  }

  inline void
  Parameter :: write ( const std :: string &filename )
  {
    // only write one parameter log file
    if( MPIManager :: rank() != 0 )
      return;

    std :: string fullname( commonOutputPath() );
    fullname += "/";
    fullname += filename;

    std :: ofstream file( fullname.c_str() );
    if( !file.is_open() )
    {
      std :: cerr << "Warning: Unable to write parameter file '"
                  << filename << "'" << std :: endl;
    }
    write( file );
    file.close();
  }

  inline void
  Parameter :: write ( std :: ostream &out )
  {
    typedef std :: map< std :: string, std :: string > :: iterator iterator;
    const iterator end = instance().params_.end();
    for( iterator it = instance().params_.begin(); it != end; ++it )
      out << it->first << ": " << it->second << std :: endl;
  }

}

#endif
