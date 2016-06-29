#include <config.h>

#include <iostream>
#include <string>
#include <tuple>

#include <dune/common/exceptions.hh>

#include <dune/fem/common/tupleforeach.hh>

int main()
try
{
  std::tuple<int,double,std::string> t{1,-0.5,"tuple"};
  Dune::Fem::for_each(t,[](const auto& entry, auto I){std::cout<<"Entry "<<I<<" : "<<entry<<std::endl;});
  return 0;
}
catch(const Dune::Exception& exception)
{
  std::cerr<<exception<<std::endl;
  return 1;
}
