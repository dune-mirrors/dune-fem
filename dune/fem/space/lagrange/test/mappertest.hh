#ifndef DUNE_FEM_LAGRANGEMAPPER_TEST_HH
#define DUNE_FEM_LAGRANGEMAPPER_TEST_HH

#include <string>

#include <dune/fem/gridpart/leafgridpart.hh>

namespace Dune
{

  namespace Fem
  {

    template< class Grid >
    class LagrangeMapper_Test
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

  } // namespace Fem

} //namespace Dune

#endif // DUNE_FEM_LAGRANGEMAPPER_TEST_HH
