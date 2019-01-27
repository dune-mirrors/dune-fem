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
    // works: auto param = Dune::Fem::parameterMap({{"test","test"},{"one","1"},{"one.one","1.1"}});
#if 0
    // fails because initializer lists can't be used with variadic templates?
    auto param = Dune::Fem::parameterDict("",{"test","test"},{"one",1},{"one.one",1.1});
#elif 0
    // first version in parameter.hh requires make_pair
    // would like to avoid the make pair
    auto param = Dune::Fem::parameterDict("",
        std::make_pair("test","test"),
        std::make_pair("one",1),
        std::make_pair("one.one",1.1));
#else
    auto param = Dune::Fem::parameterDict("","test","test","one",1,"one.one",1.1);
#endif
    std::cout << param.getValue<std::string>("test") << " "
              << param.getValue<int>("one") << " "
              << param.getValue<double>("one.one") << std::endl;
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
