#ifndef DUNE_LAGRANGEMAPPER_TEST_HH
#define DUNE_LAGRANGEMAPPER_TEST_HH

#include <string>

#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/misc/test.hh>

namespace Dune
{

  template< class Grid >
  class LagrangeMapper_Test
  : public Test
  {
  public:
    typedef Grid GridType;

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
