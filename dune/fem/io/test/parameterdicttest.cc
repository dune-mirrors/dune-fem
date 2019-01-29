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

// main programm
int main(int argc, char** argv)
{
  try{
    double a = 2;
    auto param = Dune::Fem::parameterDict("prefix.",
                    "test","test",  "one",1,  "one.one",1.1,
                    "lambdaS",[&a](){ return "a="+std::to_string(int(a)); },
                    "lambdaD",[&a](){ return a; }
                    );
    bool pass = true;
    pass &= ( param.exists("prefix.test") &&
              param.exists("prefix.one") &&
              param.exists("prefix.one.one") &&
              param.exists("prefix.lambdaS") &&
              param.exists("prefix.lambdaD") &&
             !param.exists("prefix.notfound") &&
             !param.exists("notfound") );
    pass &= param.getValue<std::string>("prefix.test") == "test";
    pass &= param.getValue<int>("prefix.one") == 1;
    pass &= param.getValue<double>("prefix.one.one") == 1.1;
    pass &= param.getValue<std::string>("prefix.lambdaS") == "a=2";
    pass &= param.getValue<double>("prefix.lambdaD") == 2;
    try{
      std::cout << param.getValue<std::string>("prefix.notfound");
      pass = false;
    } catch (Dune::Fem::ParameterNotFound &e) { }
    try{
      std::cout << param.getValue<std::string>("notfound");
      pass = false;
    } catch (Dune::Fem::ParameterNotFound &e) { }

    pass &= param.getValue<std::string>("prefix.test","testing") == "test";
    pass &= param.getValue<int>("prefix.one",2) == 1;
    pass &= param.getValue<double>("prefix.one.one",2.2) == 1.1;
    pass &= param.getValue<std::string>("prefix.lambdaS","a=5.") == "a=2";
    pass &= param.getValue<double>("prefix.lambdaD",5.) == 2;
    pass &= param.getValue<std::string>("prefix.notfound","not_found") == "not_found";
    pass &= param.getValue<std::string>("notfound","not_found") == "not_found";

    a = 4;
    pass &= param.getValue<std::string>("prefix.lambdaS") == "a=4";
    pass &= param.getValue<double>("prefix.lambdaD") == 4;

    return pass ? 0 : 1;
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
