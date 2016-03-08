#ifndef DUNE_FEM_PARAMETER_HH
#define DUNE_FEM_PARAMETER_HH

#include <fstream>
#include <iostream>
#include <string>

#include <dune/fem/io/io.hh>
#include <dune/fem/io/parameter/exceptions.hh>
#include <dune/fem/io/parameter/container.hh>
#include <dune/fem/io/parameter/reader.hh>

namespace Dune
{

  namespace Fem
  {

    /** \addtogroup Parameter

        Handling Parameters, i.e., values that can be set after compilation, in
        dune-fem is extremely easy. Just add
        \code
        Dune::MPIManager::initialize( argc, argv );
        Dune::Parameter::append( argc, argv );
        \endcode
        at the head of your main function. Parameters are strings of the
        format "key: value". Any command line argument containing a colon
        will thus be interpreted as a parameter.

        \note All parameter keys are case sensitive.
        \note Found parameters are removed from the command line.

        There are 8 static methods in Parameter to obtain the value of a
        parameter. They are divided by the following criteria:
        - \b Default \b Value: The methods taking a default value will return the
          default, if the parameter has not been specified by the user.
          Additionally, they will add "\e key : \e value" to the database. The
          methods not taking a default value will throw an exception, if the
          parameter could not be found in the database.
        - \b Return \b Value: For convenience, there is always a method (called
          getValue) returning the value of the parameter. If you do not want to
          rely on return value optimization, use the method (called get) taking
          a reference to the return variable as an argument.
        - \b Validation: It is often necessary to make sure the value satisfies
          certain constraints. For this purpose, some methods take a validator
          as an argument.
        .

        Of course, you don't have to pass every parameter on the command line.
        They can also be gathered in files. Parameter provides a kind of include
        mechanism. Whenever a parameter with key "paramfile" is encountered, the
        value is interpreted as a paramter file to include.

        Fem parameters can also be appended to the parameter list via an DGF file.
        The method
        \code
        appendDGF( "macrogrid.dgf" );
        \endcode
        reads in the DGF-block:
        \code
        FemParameter

        #
        \endcode
        from file "macrogrid.dgf". All parameters defined within this block are appended
        to the parameter list.


        If a parameter is defined multiply, the first definition is added to the
        database. All later definitions are ignored. Therefore it is important to
        know the exact behaviour of the "paramfile" parameter:
        - All parameters in the current file (or the command line) are added
          first.
        - If there were includes, they are processed in the order of appearance.
        - Should an included file have includes, they are added depth-first,
          i.e., The includes of one included file are parsed down to the last
          file included, before the includes of the next included file are
          considered.
        .

        All parameter names defined by dune-fem should conform to the following
        naming convention:
        \code
        fem.<group>.<parameter>
        \endcode
        The group name can be omitted if necessary.

        An example is the parameter
        \code
        fem.verboserank
        \endcode
        This can beused throughout the the program; by calling:
        \code
        Parameter::verbose()
        \endcode
        If verbose is set, information concerning the parameters read will
        be output to stdout.

        \b Parameter Substitution: \b
        Parameter can consist of parts of other parameters. In order to resolve this
        dependency, the  substituted parameter is put into $(NewParameter) brackets.
        A smale example for this is
        \code
        N: 128
        parameter1: macrogrid_$(N).dgf
        \endcode
        results in
        \code
        parameter1: macrogrid_128.dgf
        \endcode
        This can be used when on parameter controls several other parameters, such like the
        number of cells in the macrogrid.

        \b Shell Program executions:\b
        Out of the Parameter file shell scripts/commands can be called in order to
        calculate the value of a parameter. The object which should be executed
        has to be put into $[command] brackets.
        \code
        parameter1: $[ ./script.sh]
        \endcode
        with the script.sh
        \code
        #!/bin/bash
        echo 'HalloWorld';
        \endcode
        This smale example resolves the parameter to have the value 'HalloWorld'

        If $ is used explicite in a parameter value, $$ kills the substitution
        of the parameter.

        Here is an example usage:
        \code
        #include <dune/fem/io/parameter.hh>

        int globalFlag;

        int main ( int argc, char **argv )
        {
          Dune::Fem::MPIManager::initialize( argc, argv );
          Dune::Fem::Parameter::append( argc, argv );

          // get parameters
          double startTime = Parameter::getValue< double >( "starttime", 0.0 );
          auto validator = [ startTime ]( double time ) { return time > startTime; };
          double endTime = Parameter::getValidValue< double >( "endtime", validator );

          Dune::Fem::Parameter::get( "flag", 0, globalFlag );

          if( Dune::Fem::Parameter::verbose() )
            std::cout << "Computing from " << startTime << " to " << endTime << std::endl;

          // ...

          std::ofstream results( Dune::Fem::Parameter::outputPrefix() + "/results" );

          // ...

          Dune::Fem::Parameter::write( "parameter.log" );
        }
        \endcode

        \b Parameter deprecation \b
        Sometimes parameter names are changed and it can be difficult to find all occurrences of the old parameter
        in the code. To help with the transition it is possible to deprecate parameters by adding the line
        \code
        deprecated: old_name
        \endcode
        This will cause an exception to be thrown by any attempt in the code to access the value of the parameter
        using the old name.
     */



    // Parameter
    // ---------

    /** \class Parameter
        \brief Container for User Specified Parameters

        The class Parameter provides parameters collected from given parameter
        files or the command line in a unified and easy to use way.

        This class adheres to the singleton concept, i.e., all methods are static
        and internally use a single instance of this object to store all data.
     */
    class Parameter
    {
    public:
      static ParameterContainer &container ()
      {
        static ParameterContainer container;
        return container;
      }

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
      static void append ( int &argc, char **argv ) { container().append( argc, argv ); }

      /** \brief add a single parameter to the container
       *
       *  \param[in]  key    key of the parameter to add
       *  \param[in]  value  value of the parameter to add
       */
      static void append ( const std::string &key, const std::string &value ) { container().append( key, value ); }

      /** \brief add parameters from a file to the container
       *
       * \param[in]  filename  name of the file containing the parameters
       */
      static void append ( const std::string &filename ) { container().append( filename ); }

      /** \brief add parameters from a DGF file to the container
       *
       *  Parameters can also be read from a DGF file containing a 'FemParameter'
       *  block.
       *
       *  \param[in]  filename  name of the DGF file containing the parameters
       */
      static void appendDGF ( const std::string &filename ) { container().appendDGF( filename ); }

      /** \brief find out, whether a parameter is defined in the container
       *
       *  \param[in]   key    name of the parameter to check
       *
       *  \returns \b true, if the parameter is found in the container,
       *           \b false otherwise
       */
      static bool exists ( const std::string &key ) { return container().exists( key ); }

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
        container().get( key, value );
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
      static void get ( const std::string &key, const T &defaultValue, T &value )
      {
        container().get( key, defaultValue, value );
      }

      /** \brief get an optional parameter from the container special case for string
       *
       *  \note This method returns a default value, if the parameter cannot be
       *        found.
       *
       *  \param[in]   key           name of the parameter to get
       *  \param[in]   defaultValue  default value for this parameter
       *  \param[out]  value         value of the parameter
       */
      static void get ( const std::string &key, const char* defaultValue, std::string &value )
      {
        container().get( key, defaultValue, value );
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
      static void getValid ( const std::string &key, const Validator &validator, T &value )
      {
        container().getValid( key, validator, value );
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
      static void getValid ( const std::string &key, const T &defaultValue, const Validator &validator, T &value )
      {
        container().getValid( key, validator, value );
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
      static T getValue ( const std::string &key )
      {
        return container().getValue< T >( key );
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
        return container().getValue( key, defaultValue );
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
        return container().getValidValue( key, validator );
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
        return container().getValidValue( key, defaultValue, validator );
      }

      template< int n >
      static int getEnum ( const std::string &key, const std::string (&values)[ n ] )
      {
        return container().getEnum( key, values );
      }

      template< int n >
      static int getEnum ( const std::string &key, const std::string (&values)[ n ], int defaultValue )
      {
        return container().getEnum( key, values, defaultValue );
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
      static std::string commonOutputPath ()
      {
        return container().commonOutputPath();
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

      /** \brief obtain common input path
       *
       *   This common input path could for example be used for the
       *   determination of the location of parameter and grid files
       *
       *  \returns value of parameter 'fem.input.prefix', which defaults to '.'.
       */
      static std::string commonInputPath ()
      {
        return container().commonInputPath();
      }

      /** \brief obtain the cached value for fem.verbose */
      static bool verbose () { return container().verbose(); }

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
       *  \param[in]  fileextension  file extension if you want to have time stemps in the log files,
       *              ! must contain the "." dot for fileextensions !
       *  \param[in]  writeAll default is true
       */
      static void write ( const std::string &filename, const std::string &fileextension ="", bool writeAll = true );

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
      static void write ( std::ostream &out, bool writeAll = true ) { return container().write( out, writeAll ); }

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
       *  \param[in]  fileextension  chosen fileextension
       *  \param[in]  writeAll default is true
       */
      static void write ( const std::string &path, const std::string &filename, const std::string &fileextension, bool writeAll = true );
    };



    // Implementation of Parameter
    // ---------------------------

    inline std::string Parameter::outputPath ()
    {
      std::string path = commonOutputPath();
      if( path.empty() )
        path = "./";
      if( path[ path.length()-1 ] != '/' )
        path += '/';
      return path + "p" + std::to_string( MPIManager::rank() );
    }


    inline void Parameter::write ( const std::string &filename, const std::string &fileextension, bool writeAll )
    {
      // only write one parameter log file
      // to the common path
      if( MPIManager::rank() == 0 )
        write( commonOutputPath(), filename, fileextension, writeAll );
    }


    inline void Parameter::write ( const std::string &path, const std::string &filename, const std::string &fileextension, bool writeAll )
    {
      //create path if it does not exist
      if( !directoryExists( path ) )
        createDirectory( path );

      std::string fullname( path );
      fullname += "/";
      fullname += filename;
      fullname += fileextension;

      std::ofstream file( fullname );
      if( file.is_open() )
        write( file, writeAll );
      else
        std::cerr << "Warning: Unable to write parameter file '" << fullname << "'" << std::endl;
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
        return new ParamDefault( static_cast< const ParamDefault & >( *this ) );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PARAMETER_HH
