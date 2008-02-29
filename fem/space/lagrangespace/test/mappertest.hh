#ifndef DUNE_LAGRANGEMAPPER_TEST_HH
#define DUNE_LAGRANGEMAPPER_TEST_HH

#define LAGRANGE_CHECK_MAPPING
#define TEST_SECOND_ORDER
#define TEST_THIRD_ORDER


#include <string>
#include <dune/fem/misc/test.hh>

namespace Dune {

  class LagrangeMapper_Test : public Test {
  public:
    typedef LeafGridPart< GridType > GridPartType;
      
    LagrangeMapper_Test( std :: string gridFile )
      : gridFile_( gridFile )
    {
    }

    virtual void run();

  private:
    template< class SpaceType >
    void checkDiscreteFunction( const SpaceType &space );

    std :: string gridFile_;
  };
    
}

#endif // DUNE_LAGRANGEMAPPER_TEST_HH
