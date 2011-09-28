#ifndef DUNE_FEM_PARAMETER_HH
#define DUNE_FEM_PARAMETER_HH

#include <map>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>

#include <dune/common/exceptions.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/io/io.hh>
#include <dune/fem/misc/validator.hh>
#include <dune/fem/misc/mpimanager.hh>

namespace Dune
{
  
  /** \addtogroup Parameter
   *
   *  Handling Parameters, i.e., values that can be set after compilation, in
   *  dune-fem is extremely easy. Just add
   *  \code
   *  Dune::MPIManager::initialize( argc, argv );
   *  Dune::Parameter::append( argc, argv );
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
   *  Parameter::verbose()
   *  \endcode
   *  If verbose is set, information concerning the parameters read will
   *  be output to stdout.
   *
   *  \b Parameter Substitution: \b
   *  Parameter can consist of parts of other parameters. In order to resolve this
   *  dependency, the  substituted parameter is put into $(NewParameter) brackets.
   *  A smale example for this is
   *  \code
   *  N: 128
   *  parameter1: macrogrid_$(N).dgf
   *  \endcode
   *  results in 
   *  \code
   *  parameter1: macrogrid_128.dgf
   *  \endcode
   *  This can be used when on parameter controls several other parameters, such like the
   *  number of cells in the macrogrid.
   *
   *  \b Shell Program executions:\b
   *  Out of the Parameter file shell scripts/commands can be called in order to
   *  calculate the value of a parameter. The object which should be executed 
   *  has to be put into $[command] brackets.
   *  \code
   *  parameter1: $[ ./script.sh]
   *  \endcode
   *  with the script.sh
   *  \code
   *  #!/bin/bash
   *  echo 'HalloWorld';
   *  \endcode
   *  This smale example resolves the parameter to have the value 'HalloWorld'
   *
   *  If $ is used explicite in a parameter value, $$ kills the substitution
   *  of the parameter.
   *
   *  Here is an example usage:
   *  \code
   *  #include <dune/fem/io/parameter.hh>
   *  
   *  using Dune::Parameter;
   *
   *  int globalFlag;
   *
   *  int main ( int argc, char **argv )
   *  {
   *    Dune::MPIManager::initialize( argc, argv );
   *    Parameter::append( argc, argv );
   *
   *    // get parameters
   *    double startTime = Parameter::getValue< double >( "starttime", 0.0 );
   *    Dune::ValidateGreater< double > validator( startTime );
   *    double endTime = Parameter::getValidValue< double >( "endtime", validator );
   *    Parameter::get( "flag", 0, globalFlag );
   *    
   *    if( Parameter::verbose() )
   *      std::cout << "Computing from " << startTime << " to " << endTime << std::endl;
   *    
   *    // ...
   *
   *    std::ofstream results( (Parameter::outputPrefix() + "/results").c_str() );
   *
   *    // ...
   *
   *    Parameter::write( "parameter.log" );
   *  }
   *  \endcode
   */

  class ParameterNotFound
  : public Exception
  {};

  class ParameterInvalid
  : public Exception
  {};



  // ParameterParser
  // ---------------

  template< class T >
  struct ParameterParser
  {
    static bool parse ( const std::string &s, T &value )
    {
      std::istringstream in( s );
      in >> value;
      if( in.fail() )
        return false;
      char eof;
      in >> eof;
      return in.eof();
    }

    static std::string toString ( const T &value )
    {
      std::ostringstream out;
      out << value;
      return out.str();
    }
  };

  template<>
  struct ParameterParser< bool >
  {
    static bool parse ( const std::string &s, bool &value )
    {
      std::string w;
      if( ParameterParser< std::string >::parse( s, w ) )
      {
        if( (w == std::string( "false" )) || (w == std::string( "no" )) || (w == std::string( "0" )) )
        {
          value = false;
          return true;
        }

        if( (w == std::string( "true" )) || (w == std::string( "yes" )) || (w == std::string( "1" )) )
        {
          value = true;
          return true;
        }
      }
      return false;
    }

    static std::string toString ( const bool &value )
    {
      return std::string( value ? "true" : "false" );
    }
  };

  template< class F, int m, int n >
  struct ParameterParser< FieldMatrix< F, m, n > >
  {
    static bool parse ( const std::string &s, FieldMatrix< F, m, n > &value )
    {
      std::istringstream in( s );
      char c;
      for( int i = 0; i < m; ++i )
      {
        if( i > 0 )
        {
          in >> c;
          if( c != ',' )
            return false;
        }

        for( int j = 0; j < n; ++j )
          in >> value[ i ][ j ];
      }
      in >> c; // read eof
      return in.eof();
    }

    static std::string toString ( const FieldMatrix< F, m, n > &value )
    {
      std::ostringstream out;
      for( int i = 0; i < m; ++i )
      {
        out << (i > 0 ? "," : "");
        for( int j = 0; j< n; ++j )
          out << " " << value[ i ][ j ];
      }
      return out.str();
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

    struct Value
    {
      enum ShadowStatus{ unresolved, resolved, resolving};
      Value() : used(false), hasDefault(false), isDefault(false), shadowStatus(unresolved) {}

      std::string value, fileName, defaultValue;
      bool used, hasDefault, isDefault;
      ShadowStatus shadowStatus;
    };

    struct DGFBlock;

    enum CheckDefaultType { checkDefaultDisable, checkDefaultEnable };

    typedef std::map< std::string, Value > ParameterMapType;

    Parameter ()
    : verboseRank_( -1 ),
      enableShadows_( false )
    {}

    // Prohibit copying and assignment
    Parameter ( const ThisType & );
    ThisType &operator= ( const ThisType & );

    static ThisType &instance ()
    {
      static ThisType theInstance;
      return theInstance;
    }

    Value *find ( const std::string &key );

    const std::string &map ( const std::string &key, const CheckDefaultType checkDefault = checkDefaultEnable );
    template <class T>
    const std::string &map ( const std::string &key, const T &value );

    const std::string &insert ( const std::string &key, const std::string &value );

    static std::string stripComment ( const std::string &line );
    static std::string trim ( const std::string &s );

    bool insert ( const std::string &s, std::queue< std::string > &includes );

    void processDGF ( const std::string &filename );
    void processFile ( const std::string &filename );
    void processIncludes( std::queue< std::string > &includes );

    std::string resolveEscape ( const std::string &key, std::string &value );
    void resolveShadows( const std::string &key, Value &val );
    std::string getShadowKey( const std::string key, const char delimter, std::string &value );

#if 0
    void replace ( const std::string &key, const std::string &value );
#endif

  public:
    /** \brief add parameters from the command line
        RangeType gRight;
     *
     *  This mehtod adds all parameters (strings containing a colon) in the
     *  command line to the container. The parameters are then removed from the
     *  command line.
     * 
     * \param[in]  argc  number of arguments (as given to main)
     * \param[in]  argv  vector of arguments (as given to main)
     */
    static void append ( int &argc, char **argv );

    /** \brief add a single parameter to the container
     *
     *  \param[in]  key    key of the parameter to add
     *  \param[in]  value  value of the parameter to add
     */
    static void append ( const std::string &key, const std::string &value )
    {
      if( key != "paramfile" )
      {
        instance().curFileName_ = "program code";
        instance().insert( key, value );
      }
      else
        instance().processFile( value );
    }

    /** \brief add parameters from a file to the container
     * 
     * \param[in]  filename  name of the file containing the parameters
     */
    static void append ( const std::string &filename )
    {
      instance().processFile( filename );
    }

    /** \brief add parameters from a DGF file to the container
     *
     *  Parameters can also be read from a DGF file containing a 'FemParameter'
     *  block.
     * 
     *  \param[in]  filename  name of the DGF file containing the parameters
     */
    static void appendDGF ( const std::string &filename )
    {
      instance().processDGF( filename );
    }

    /** \brief clear all parameters
     */
    static void clear ()
    {
      instance().params_.clear();
    }
   
#if 0
    template <class T>
    static void replaceKey ( const std::string& key, const T& value )
    {
      std::stringstream valueStr;
      valueStr << value;
      instance().replace( key,valueStr.str() );
    }
#endif

    /** \brief find out, whether a parameter is defined in the container
     *
     *  \param[in]   key    name of the parameter to check
     *
     *  \returns \b true, if the parameter is found in the container,
     *           \b false otherwise
     */
    static bool exists ( const std::string &key )
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
    static void get ( const std::string &key, T &value )
    {
      if( !ParameterParser< T >::parse( instance().map( key ), value ) )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
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
    static void get ( const std::string &key, const T &defaultValue, T &value );

    /** \brief get an optional parameter from the container special case for string
     *
     *  \note This method returns a default value, if the parameter cannot be
     *        found.
     *
     *  \param[in]   key           name of the parameter to get
     *  \param[in]   defaultValue  default value for this parameter
     *  \param[out]  value         value of the parameter
     */
    static void get ( const std::string &key, const char* defaultValue, std::string &value );

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
    static void getValid ( const std::string &key, const Validator &validator, T &value );

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
    static void getValid ( const std::string &key, const T &defaultValue, const Validator &validator, T &value );
    
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
    static T getValue ( const std::string &key )
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
    static T getValue ( const std::string &key, const T &defaultValue )
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
    static T getValidValue ( const std::string &key, const Validator &validator )
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
    static T getValidValue ( const std::string &key, const T &defaultValue, const Validator &validator )
    {
      T value;
      getValid( key, defaultValue, validator, value );
      return value;
    }

    template< int n >
    static int getEnum ( const std::string &key, const std::string (&values)[ n ] );
    template< int n >
    static int getEnum ( const std::string &key, const std::string (&values)[ n ], const int defaultValue );
  protected:
    template< int n >
    static int getEnumeration( const std::string &key, const std::string& value, const std::string (&values)[ n ]);

  public:  
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
    static std::string commonOutputPath ()
    {
      return instance().map( "fem.prefix", std::string(".") );
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
    static std::string outputPath ();

    /** \brief obtain the value for fem.prefix defaults to '.'
     */
    DUNE_DEPRECATED static const std::string &prefix ()
    {
      return instance().map( "fem.prefix", std::string(".") );
    }

    /** \brief obtain the cached value for fem.verbose
     */
    static bool verbose ()
    {
      return (instance().verboseRank_ == MPIManager::rank());
    }
  
    /** \brief write the parameter database to a file
     *
     *  This method writes paramters to the given file.
     *  If the second parameter is true all parameters are written;
     *  otherwise only used parameters which do not coincide with the default value
     *  are written.
     *
     *  \note This method is safe for parallel jobs. Parameters are only written
     *        on rank 0.
     *
     *  \param[in]  filename  name of the file to store the parameters in; prefix() is used.
     *  \param[in]  writeAll default is true
     */
    static void write ( const std::string &filename, bool writeAll = true );
    /** \brief write the parameter database to a stream
     *
     *  This method writes paramters to the given stream.
     *  If the second parameter is true all parameters are written;
     *  otherwise only used parameters which do not coincide with the default value
     *  are written.
     *
     *  \note This method is \b not safe for parallel jobs. Parameters are
     *        written on all ranks. 
     *
     *  \param[in]  out stream for the parameters.
     *  \param[in]  writeAll default is true
     */
    static void write ( std::ostream &out, bool writeAll = true );

  protected:  
    friend class PersistenceManager ;

    /** \brief write the parameter database to a file
     *
     *  This method writes paramters to the given file.
     *  If the second parameter is true all parameters are written;
     *  otherwise only used parameters which do not coincide with the default value
     *  are written.
     *
     *  \note This method is safe for parallel jobs. Parameters are only written
     *        on rank 0.
     *
     *  \param[in]  path  path where filename is stored (for parallel writing)
     *  \param[in]  filename  name of the file to store the parameters in; prefix() is used.
     *  \param[in]  writeAll default is true
     */
    static void write ( const std::string& path, const std::string &filename, 
                        bool writeAll = true  );
  private:
    std::string curFileName_;
    int curLineNumber_;
    ParameterMapType params_;
    int verboseRank_;
    bool enableShadows_;
  };



  // Parameter::DGFBlock
  // -------------------

  struct Parameter::DGFBlock
  : dgf::BasicBlock
  {
    DGFBlock ( std::istream &in )
    : BasicBlock( in, "FemParameter" )
    {}

    bool advance () { return getnextline(); }
    std::string getLine () const { return line.str(); }
  };



  // Private Methods
  // ---------------

  inline Parameter::Value *Parameter::find ( const std::string &key )
  {
    ParameterMapType::iterator it = params_.find( key );
    return (it != params_.end()) ? &(it->second) : 0;
  }

  inline const std::string &Parameter::map ( const std::string &key, const CheckDefaultType checkDefault )
  {
    Value *val = find( key );

    if (!val)
      DUNE_THROW( ParameterNotFound, "Parameter '" << key << "' not found." ); 

    resolveShadows( key, *val );


    if( checkDefault == checkDefaultEnable )
    {
      if( val->used && val->hasDefault )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' used first with and then without default." );
      val->used = true;
      val->hasDefault = false;
    }

    return val->value;
  }

  template< class T >
  inline const std::string &
  Parameter::map ( const std::string &key, const T &value )
  {
    Value insVal;
    insVal.value = ParameterParser< T >::toString( value );
    insVal.fileName = "default";
    insVal.used = false;
    insVal.shadowStatus = (enableShadows_ ? Value::unresolved : Value::resolved);

    std::pair< ParameterMapType::iterator, bool > info
      = params_.insert( std::make_pair( key, insVal ) );
    Value &val = info.first->second;

    if ( info.second )  
    {
      if ( verbose() )
        std::cout << "Adding default: " << key << " = " << value << std::endl;
      val.isDefault = true;
    }
    else if ( val.used )
    {
      if (!val.hasDefault ) 
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' used first without and then with default." );
      if ( val.defaultValue != insVal.value )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' used with different default values." );
    }
    else
    {
      T setValue( value );
      if( ParameterParser< T >::parse( val.value, setValue ) )
        val.isDefault = (setValue == value);
    }
    val.used = true;
    val.hasDefault = true;
    val.defaultValue = insVal.value;

    resolveShadows( key, val );

    return val.value;
  }

  inline const std::string &
  Parameter::insert ( const std::string &key, const std::string &value )
  {
    Value insVal;
    insVal.value = value;
    insVal.fileName = curFileName_;
    insVal.used = false;
    insVal.shadowStatus = (enableShadows_ ? Value::unresolved : Value::resolved);

    std::pair< ParameterMapType::iterator, bool > info
      = params_.insert( std::make_pair( key, insVal ) );
    Value &val = info.first->second;

    if( verbose() )
    {
      std::cout << curFileName_ << "[" << curLineNumber_ << "]: ";
      if( info.second )
        std::cout << "Adding " << key << " = " << value << std::endl;
      else 
        std::cout << "Ignored " << key << " = " << value
                  << ", using " << val.value << std::endl;
    }

    return val.value;
  }

  inline std::string Parameter::stripComment ( const std::string &line )
  {
    size_t size = line.size();
    size_t end = line.find_first_of ( "%#$" );

    while( (end != std::string::npos) && (line[end] =='$') )
    {
      if( end+2 < size ) 
        end = line.find_first_of ( "%#$", end+2 );
      else
        end = std::string::npos;
    }

    return trim( line.substr( 0, end ) );
  }

  inline std::string Parameter::trim( const std::string &s )
  {
    size_t first = s.find_first_not_of(" \t\n");
    size_t last = s.find_last_not_of(" \t\n");

    if( first == std::string::npos )
    {
      assert( last == std::string::npos );
      return std::string("");
    }
    else
    {
      assert( last != std::string::npos );
      return s.substr(first, last - first + 1);
    }
  }


  inline bool
  Parameter::insert ( const std::string &s, std::queue< std::string > &includes )
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

    unsigned int value_end = value_start;
    for( unsigned int i = 0; i < size; ++i )
    {
      if( (s[ i ] != ' ') && (s[ i ] != '\t') )
        value_end = i+1;
    }

    if( value_start >= size )
      return false;

    std::string key = s.substr( key_start, key_end - key_start );
    std::string value = s.substr( value_start, value_end - value_start );

    if( key != "paramfile" )
    {
      const std::string &actual_value = insert( key, value );
      if( key == "fem.verbose" )
      {
        bool verbose;
        ParameterParser< bool >::parse( actual_value, verbose );
        verboseRank_ = (verbose ? MPIManager::rank() : -1);
        std::cout << "Warning: Using deprecated parameter 'fem.verbose'; "
                    << "use 'fem.verboserank' instead." << std::endl;
      }
      if( key == "fem.verboserank" )
      {
        ParameterParser< int >::parse( actual_value, verboseRank_ );
        if( (verboseRank_ < -1) || (verboseRank_ >= MPIManager::size() ) )
        {
          std::cout << "Warning: Parameter 'fem.verboserank' is neither a "
                      << "valid rank nor -1." << std::endl;
        }
      }
      if( key == "fem.resolvevariables" )
         ParameterParser< bool >::parse( actual_value, enableShadows_ );
      }
    else
      includes.push( value );
    return true;
  }


  inline void
  Parameter::processDGF ( const std::string &filename )
  {
    if( verbose() )
      std::cout << "Parameter: Processing DGF '" << filename << "'..."
                   << std::endl;

    std::ifstream file( filename.c_str() );
    if( file.is_open() )
    {
      curFileName_ = filename;
      curLineNumber_ = 0;
      std::queue< std::string > includes;
      if( DuneGridFormatParser::isDuneGridFormat( file ) )
      {
        DGFBlock block( file );
        if( block.isactive() )
        {
          while( block.advance() )
          {
            ++curLineNumber_;
            const std::string line = stripComment( block.getLine() );
            if( line.size() == 0 )
              continue;
            insert( line, includes );
          }
        }
      }
      processIncludes( includes );
    }
    else
      std::cerr << "Warning: Unable to read DGF file '" << filename << "'" << std::endl;
  }


  inline void
  Parameter::processFile ( const std::string &filename )
  {
    if( verbose() )
      std::cout << "Parameter: Processing '" << filename << "'..."
                   << std::endl;
    
    std::ifstream file( filename.c_str() );
    if( file.is_open() )
    {
      curFileName_ = filename;
      curLineNumber_ = 0;
      std::queue< std::string > includes;
      while( !file.eof() )
      {
        std::string line;
        std::getline( file, line );
        curLineNumber_++;
        line = stripComment( line );
        if( line.size() == 0 )
          continue;
        insert( line, includes );
      }
      file.close();

      processIncludes( includes );
    }
    else
      std::cerr << "Warning: Unable to read parameter file '"
                  << filename << "'" << std::endl;
  }

  inline void
  Parameter::processIncludes( std::queue< std::string > &includes )
  {
    while( !includes.empty() )
    {
      Value val;
      val.value = includes.front();
      includes.pop();

      val.shadowStatus = ( enableShadows_ ? Value::unresolved : Value::resolved );
      resolveShadows( "paramfile", val );

      processFile( val.value );
    }
  }

#if 0
  inline void
  Parameter::replace ( const std::string &key, const std::string &value )
  {
    if( verbose() )
      std::cout << "Parameter: Replacing " << key << " = " << value
                  << std::endl;
    params_[ key ] = value;
  }
#endif


  // Public methods
  // --------------

  inline void
  Parameter::append ( int &argc, char **argv )
  {
    std::queue< std::string > includes;
    instance().curFileName_ = "programm arguments";
    int &i = instance().curLineNumber_;
    for( i = 1 ; i < argc; ++i )
    {
      if( !instance().insert( std::string( argv[ i ] ), includes ) )
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
  Parameter::get ( const std::string &key, const T &defaultValue, T &value )
  {
    if( !ParameterParser< T >::parse( instance().map( key, defaultValue ), value ) )
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
  }

  inline void
  Parameter::get ( const std::string &key, const char* defaultValue, std::string &value )
  {
    get(key,std::string(defaultValue),value);
  }

  template< class T, class Validator >
  inline void
  Parameter::getValid ( const std::string &key,
                        const Validator &validator,
                        T &value )
  {
    if( !ParameterParser< T >::parse( instance().map( key ), value ) )
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
    if( !validator( value ) )
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
  }

  template< class T, class Validator >
  inline void
  Parameter::getValid ( const std::string &key,
                        const T &defaultValue,
                        const Validator &validator,
                        T &value )
  {
    bool valid = true ;
    if( !ParameterParser< T >::parse( instance().map( key, defaultValue ), value ) ) valid = false;
    if( !validator( value ) ) valid = false;

    if( ! valid ) 
    {
      std::cerr << std::endl << "Parameter '" << key << "' invalid." << std::endl;
      validator.print( std::cerr );
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
    }
  }
  
  inline std::string
  Parameter::outputPath ()
  {
    const std::string &prefix = instance().map( "fem.prefix", std::string(".") );
    std::ostringstream path;
    if( !prefix.empty() )
    {
      path << prefix;
      if( prefix[ prefix.length()-1 ] != '/' )
        path << '/';
    }
    else
      path << "./";
    path << 'p' << MPIManager::rank();
    return path.str();
  }


  template< int n >
  inline int
  Parameter::getEnum ( const std::string &key, const std::string (&values)[ n ] )
  {
    return getEnumeration( key, instance().map( key ), values );
  }
  template< int n >
  inline int
  Parameter::getEnum ( const std::string &key, const std::string (&values)[ n ], const int defaultValue )
  {
    return getEnumeration( key, 
                           instance().map( key, values[ defaultValue ] ),
                           values);
  }

  template< int n >
  inline int
  Parameter::getEnumeration ( const std::string &key, 
                              const std::string& value, 
                              const std::string (&values)[ n ] )
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
      {
        std::cerr << values[ i ];
        if( i < n-1 ) std::cerr << ", "; 
      }
      std::cerr << std::endl << std::endl; 
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
    }
    return j;
  }


  inline std::string
  Parameter::resolveEscape ( const std::string &key, std::string &value )
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
      return map( getShadowKey( key, ')', value ), checkDefaultDisable );

    case '[':
      return trim( executeCommand( getShadowKey( key, ']', value ) ) );

    default:
      DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );
    }
  }



  inline std::string 
  Parameter::getShadowKey( const std::string key, const char delimiter, std::string &value )
  {
    std::string shadowKey; 

    while( true )
    {
      size_t startPoint = value.find_first_of( std::string("$") + delimiter );        
      
      if( startPoint == std::string::npos ) 
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid." );

      shadowKey += value.substr( 0, startPoint );
      const char startChar = value[ startPoint ];

      value.replace(0, startPoint+1, "" );

      if( startChar == delimiter ) 
        return shadowKey;
      assert( startChar == '$' );

      shadowKey += resolveEscape( key, value );
    }
  }

  inline void 
  Parameter::resolveShadows ( const std::string &key, Value &val )
  {
    std::string &realValue = val.value;
    if( val.shadowStatus != Value::resolved )
    {      
      if ( val.shadowStatus == Value::resolving )
        DUNE_THROW( ParameterInvalid, "Parameter '" << key << "' invalid, contains infinite loop." );

      val.shadowStatus = Value::resolving;
      std::string realValueHelper;
      realValue.swap(realValueHelper);

      while( !realValueHelper.empty() )
      {
        size_t startPoint = realValueHelper.find_first_of('$');        
        realValue += realValueHelper.substr( 0, startPoint );

        if( startPoint == std::string::npos ) 
          break;

        realValueHelper.replace( 0, startPoint+1, "" );

        realValue += resolveEscape( key, realValueHelper );
      }
      val.shadowStatus = Value::resolved;
    }
  }

  inline void
  Parameter::write ( const std::string &filename, bool writeAll )
  {
    // only write one parameter log file
    // to the common path 
    if( MPIManager::rank() != 0 ) return;

    write( commonOutputPath(), filename, writeAll );
  }

  inline void
  Parameter::write ( const std::string &path, 
                     const std::string &filename , bool writeAll ) 
  {
    std::string fullname( path );
    fullname += "/";
    fullname += filename;

    std::ofstream file( fullname.c_str() );
    if( !file.is_open() )
    {
      std::cerr << "Warning: Unable to write parameter file '"
                << fullname << "'" << std::endl;
      return ;
    }

    write( file, writeAll );
    file.close();
  }

  inline void
  Parameter::write ( std::ostream &out, bool writeAll )
  {
    typedef std::map< std::string,std::map<std::string,std::string> > WriteMap;
    typedef ParameterMapType::iterator Iterator;
    typedef WriteMap::iterator WriteIterator;
    typedef std::map<std::string,std::string>::iterator ParamIterator;

    WriteMap writeMap;
    const Iterator end = instance().params_.end();
    for( Iterator it = instance().params_.begin(); it != end; ++it )
    {
      const Value &val = it->second;
      if ( writeAll || (val.used && !val.isDefault) )
        writeMap[val.fileName][ (val.used?"":"# ") + it->first] = val.value;
    }
    const WriteIterator wEnd = writeMap.end();
    for( WriteIterator wit = writeMap.begin(); wit != wEnd; ++wit )
    {
      out << "# from " << wit->first << std::endl;
      const ParamIterator end = wit->second.end();
      for( ParamIterator it = wit->second.begin(); it != end; ++it )
        out << it->first << ": " << it->second << std::endl;
      out << std::endl;
    }
  }



  // LocalParameter
  // --------------

  // Helper class for Parameter structures for classes
  template< class ParamDefault, class ParamImpl >
  struct LocalParameter
  : public ParamDefault
  {
    virtual ~LocalParameter ()
    {}

    virtual ParamDefault *clone () const
    {
      return new ParamImpl( asImp() );
    }

  protected:
    const ParamImpl &asImp () const
    {
      return static_cast< const ParamImpl & >( *this );
    }
  };

  template< class ParamDefault >
  struct LocalParameter< ParamDefault, ParamDefault >
  {
    virtual ~LocalParameter ()
    {}

    virtual ParamDefault *clone () const
    {
      return new ParamDefault();
    }
  };

}

#endif // #ifndef DUNE_FEM_PARAMETER_HH
