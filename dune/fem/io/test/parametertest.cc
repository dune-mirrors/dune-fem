#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <string>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/fvector.hh>  // definition of field vectors

#include <dune/fem/io/parameter.hh> // include parameters


// main programm
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::Fem::MPIManager :: initialize( argc, argv );

    Dune::Fem::Parameter::append(argc,argv);
    if (argc == 2) {
      Dune::Fem::Parameter::append(argv[1]);
    }
    else
    {
      Dune::Fem::Parameter::append("testparameterfile");
    }

    // get user id and date
    std::string userId, date;
    Dune::Fem::Parameter::get("user", userId );
    Dune::Fem::Parameter::get("date", date );

    // get Project discription
    std::string project;
    Dune::Fem::Parameter::get( "project", project );

    // get the parameter for the diffusion
    double diffusion = Dune::Fem::Parameter::getValidValue( "diffusion", 1.0, []( const double &val ) { return val > 0; } );

    // get the macro grid filename, using a lambda to verify 'no white spaces' in the filename
    std::string macrogridname;
    auto valid = []( const std::string &name ){ return (name.find_first_of( " \t" ) == std::string::npos); };
    Dune::Fem::Parameter::getValid("macrogrid", valid, macrogridname );

    // get velocity for the advection part
    typedef Dune::FieldVector< double, 2>  VectorType;
    VectorType velocity;
    Dune::Fem::Parameter::get("velocity", velocity );

    {
      // do some calculations
    }

    // get output path
    std::string outputPath;
    Dune::Fem::Parameter::get("outputpath", outputPath );

    // just a small example that mmultiplications are possible using shell scripts
    double multi = Dune::Fem::Parameter::getValue<double>("multiplication", 1.0 );
    double factor1 = Dune::Fem::Parameter::getValue<double>("factor1", 1.0 );
    double factor2 = Dune::Fem::Parameter::getValue<double>("factor2", 1.0 );

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

    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
