#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <iostream>
#include <fstream>
#include <string>

#include <dune/common/exceptions.hh> // We use exceptions
#include <dune/common/fvector.hh>  // definition of field vectors

#include <dune/fem/io/parameter.hh> // include parameters


using namespace Dune;
using namespace Fem;

// main programm
int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    MPIManager :: initialize( argc, argv );

    Parameter::append(argc,argv);
    if (argc == 2) {
      Parameter::append(argv[1]);
    }
    else
    {
      Parameter::append("testparameterfile");
    }

    // get user id and date
    std::string userId, date;
    Parameter::get("user", userId );
    Parameter::get("date", date );

    // get Project discription
    std::string project;
    Parameter::get( "project", project );

    // get the parameter for the diffusion
    double diffusion = Parameter::getValue<double>("diffusion", 1.0 );


    // get the macro grid filename, using the NoWhiteSpaceValidator, so that no white
    // space is in the filename
    std::string macrogridname;
    Parameter::getValid("macrogrid",  NoWhiteSpaceValidator(), macrogridname );

    // get velocity for the advection part
    typedef Dune::FieldVector< double, 2>  VectorType;
    VectorType velocity;
    Parameter::get("velocity", velocity );


    {
      // do some calculations
    }

    // get output path
    std::string outputPath;
    Parameter::get("outputpath", outputPath );

    // just a small example that mmultiplications are possible using shell scripts
    double multi = Parameter::getValue<double>("multiplication", 1.0 );
    double factor1 = Parameter::getValue<double>("factor1", 1.0 );
    double factor2 = Parameter::getValue<double>("factor2", 1.0 );

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
    Parameter::write("parameter.log");

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
