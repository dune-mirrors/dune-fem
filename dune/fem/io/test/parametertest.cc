#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // #ifdef HAVE_CONFIG_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parametertreeparser.hh>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/parameter/parametertree.hh>

// small ini data structure for testing the parameter tree reader
const char iniTree[]
  = "[test]\n"
    "string = Hello\n"
    "integer = 42\n";

// main programm
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::Fem::MPIManager::initialize( argc, argv );

    auto &parameterContainer = Dune::Fem::Parameter::container();

    parameterContainer.append(argc,argv);
    if (argc == 2) {
      parameterContainer.append(argv[1]);
    }
    else
    {
      parameterContainer.append( "testparameterfile" );
    }

    Dune::Fem::ParameterReader parameter = Dune::Fem::Parameter::container();

    // get user id and date
    std::string userId, date;
    parameter.get("user", userId );
    parameter.get("date", date );

    // get Project discription
    std::string project;
    parameter.get( "project", project );

    assert( !parameter.exists( "KEY_THAT_SHOULD_NEVER_BE_SET" ) );

    int i = -42;
    parameter.get( "KEY_THAT_SHOULD_NEVER_BE_SET", 42, i );
    if( i != 42 )
      DUNE_THROW( Dune::Exception, "Default was ignored");

    // get the parameter for the diffusion
    double diffusion = parameter.getValidValue( "diffusion", 1.0, []( const double &val ) { return val > 0; } );

    try
    {
      parameter.getValidValue( "diffusion", 1.0, []( const double &val ) { return false; } );
      DUNE_THROW( Dune::Exception, "Invalid Validator application" );
    }
    catch ( const Dune::Fem::ParameterInvalid& ) {}

    // get the macro grid filename, using a lambda to verify 'no white spaces' in the filename
    std::string macrogridname;
    auto valid = []( const std::string &name ){ return (name.find_first_of( " \t" ) == std::string::npos); };
    parameter.getValid("macrogrid", valid, macrogridname );

    // get velocity for the advection part
    typedef Dune::FieldVector< double, 2>  VectorType;
    VectorType velocity;
    parameter.get("velocity", velocity );

    // get output path
    std::string outputPath;
    parameter.get("outputpath", outputPath );

    // just a small example that mmultiplications are possible using shell scripts
    double multi = parameter.getValue<double>("multiplication", 1.0 );
    double factor1 = parameter.getValue<double>("factor1", 1.0 );
    double factor2 = parameter.getValue<double>("factor2", 1.0 );

    //*****************OutPut*********************//
    std::cout<<"User: "<< userId <<" started his compution of Project: "<< project <<"\nThe acctual date is: "<< date <<std::endl;
    std::cout<<std::endl;
    std::cout<<"Parameter for this project are:"<<std::endl;
    std::cout<<"Diffusion: \t"<< diffusion << std::endl;
    std::cout<<"Velocity:  \t"<< velocity << std::endl;
    std::cout<<"Macrogrid: \t"<< macrogridname << std::endl;
    std::cout<<std::endl;
    std::cout<<"Output path is: " << outputPath << std::endl;
    std::cout<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Multiplication between: " << factor1 << " and " <<
                factor2 <<"\ncalculated with a shell script is: "
                << multi << std::endl;

    // finallay write the parameter log
    Dune::Fem::Parameter::write("parameter.log");

    // test reading from a dune-common parameter tree
    std::cout << ">>> Testing parameter reader for ParameterTree..." << std::endl;
    Dune::ParameterTree parameterTree;
    std::istringstream iniStream( iniTree );
    Dune::ParameterTreeParser::readINITree( iniStream, parameterTree );
    const auto parameterReader = Dune::Fem::parameterReader( parameterTree );
    if( parameterReader.getValue< std::string >( "test.string" ) != "Hello" )
      DUNE_THROW( Dune::Exception, "Unable to read string parameter from parameter tree" );
    if( parameterReader.getValue< int >( "test.integer" ) != 42 )
      DUNE_THROW( Dune::Exception, "Unable to read integer parameter from parameter tree" );
    if( parameterReader.getValue< int >( "undefined.integer", 123 ) != 123 )
      DUNE_THROW( Dune::Exception, "Unable to read default integer parameter from parameter tree" );

    return 0;
  }
  catch (const Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
