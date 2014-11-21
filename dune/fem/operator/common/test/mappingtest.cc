#include <config.h>

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/fem/operator/common/mapping.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/space/common/functionspace.hh>

using namespace Dune;
using namespace Fem;

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


class C : public Mapping<
          DomainFieldType,
          RangeFieldType,
          DomainType,
          RangeType >
{
  typedef SpaceType :: DomainType DomainType;
  typedef SpaceType :: RangeType  RangeType;

public:
  void operator () (const DomainType& arg, RangeType& dest) const
  {
    dest = arg;
    std::cout << "Call C :: apply \n";
  }

};


int main ()
{
  A a;
  B b;
  C c;

  MappingType m = a + b * 2.0 + a / 2.0 + b - 3.0 * a;
  DomainType arg = 1.0;
  DomainType dest;

  MappingType m1;

  m1( arg, dest );

  m1 = c;
  m1( arg, dest );


  m1 = a + b * 2.0 + c; // result should be 4
  m1( arg, dest );
  std::cout << "Result (should be 4) of m1 combined mapping is: " << dest << std::endl;

  MappingType m2( m1 );
  m1( arg, dest );
  std::cout << "Result (should be 4) of m2 combined mapping is: " << dest << std::endl;

  m(arg,dest);

  std::cout << "Result (should be 1.5) of combined mapping is: " << dest << std::endl;
  MappingType m3 = a - b * 2.0;
  m3(arg,dest);
  std::cout << "Result (should be -1.0) of combined mapping is: " << dest << std::endl;
  return 0;
}
