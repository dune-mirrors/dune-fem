#include <config.h>
#include <iostream> 

#include <dune/common/fvector.hh>
#include <dune/fem/operator/common/mapping.hh>
#include <dune/fem/space/common/functionspace.hh>

using namespace Dune;

typedef FunctionSpace<double,double,1,1> SpaceType;
typedef SpaceType :: DomainFieldType DomainFieldType; 
typedef SpaceType :: RangeFieldType RangeFieldType; 
typedef SpaceType :: DomainType DomainType; 
typedef SpaceType :: RangeType  RangeType; 

typedef Mapping<
          DomainFieldType,
          RangeFieldType,
          DomainType,
          RangeType > MappingType ;

class A : public Mapping<
          DomainFieldType, 
          RangeFieldType, 
          DomainType, 
          RangeType >
{
  typedef SpaceType :: DomainType DomainType; 
  typedef SpaceType :: RangeType  RangeType; 

public:   
  void apply (const DomainType& arg, RangeType& dest) const 
  {
    dest = arg;
    std::cout << "Call A :: apply \n";
  }

};

class B : public Mapping<
          DomainFieldType, 
          RangeFieldType, 
          DomainType, 
          RangeType >
{
  typedef SpaceType :: DomainType DomainType; 
  typedef SpaceType :: RangeType  RangeType; 

public:   
  void apply (const DomainType& arg, RangeType& dest) const 
  {
    dest = arg;
    std::cout << "Call B :: apply \n";
  }

};


int main () 
{
  A a; 
  B b;

  MappingType m = a + b * 2.0 + a / 2.0 + b - 3.0 * a;
  DomainType arg = 1.0; 
  DomainType dest;

  m(arg,dest);

  std::cout << "Result of combined mapping is: " << dest << std::endl;
  MappingType m2 = a - b * 2.0;
  m2(arg,dest);
  std::cout << "Result of combined mapping is: " << dest << std::endl;
  return 0;
}
