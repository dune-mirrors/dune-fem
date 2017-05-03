#include "mappertest.hh"

#include <dune/fem/function/common/common.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/common/gridpart.hh>
#include <dune/fem/space/lagrange.hh>

namespace Dune
{

  namespace Fem
  {

    template< class Grid >
    void LagrangeMapper_Test< Grid >::run ()
    {
      typedef GridSelector::GridType GridType;
      static const int dimworld = GridSelector::dimworld;

      GridPtr< GridType > gridPtr( gridFile_ );
      GridType &grid = *gridPtr;
      //grid.globalRefine( 2 );
      GridPartType gridPart( grid );

      typedef FunctionSpace< double, double, dimworld, dimworld > FunctionSpaceType;

      {
        std::cout << "Testing linear function." << std::endl;
        LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > space( gridPart );
        checkDiscreteFunction( space );
      }

#ifdef TEST_SECOND_ORDER
      {
        std::cout << "Testing quadratic function." << std::endl;
        LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 > space( gridPart );
        checkDiscreteFunction( space );
      }
#endif // #ifdef TEST_SECOND_ORDER

#ifdef TEST_THIRD_ORDER
      {
        std::cout << "Testing cubic function." << std::endl;
        LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 3 > space( gridPart );
        checkDiscreteFunction( space );
      }
#endif // #ifdef TEST_THIRD_ORDER
    }



    template< class Grid >
    template< class SpaceType >
    void LagrangeMapper_Test< Grid >::checkDiscreteFunction ( const SpaceType &space )
    {
      typedef AdaptiveDiscreteFunction< SpaceType > DiscreteFunctionType;
      typedef typename SpaceType::FunctionSpaceType FunctionSpaceType;
      typedef typename DiscreteFunctionType::DofType DofType;

      static const int dimworld = SpaceType::GridPartType::dimensionworld;

      std::cout << "size of space: " << space.size() << " (= " << (space.size() / dimworld) << " * " << dimworld << ")" << std::endl;

      const auto id = gridFunctionAdapter( Identity< FunctionSpaceType >(), space.gridPart() );

      DiscreteFunctionType u( "u", space );
      u.clear();

      std::cout << std::endl;
      std::cout << "Phase I: Setting each DoF of a discrete function to its global coordinate..." << std::endl;

      typename decltype( id )::LocalFunctionType idLocal( id );
      TemporaryLocalFunction< SpaceType > vLocal( space ), wLocal( space );
      for( const auto &entity : space )
      {
        const auto &interpolation = space.interpolation( entity );

        idLocal.init( entity );
        vLocal.init( entity );
        interpolation( idLocal, vLocal.localDofVector() );

        u.setLocalDofs( entity, vLocal.localDofVector() );

        wLocal.init( entity );
        interpolation( vLocal, wLocal.localDofVector() );

        if( !std::equal( vLocal.localDofVector().begin(), vLocal.localDofVector().end(), wLocal.localDofVector().begin(), [] ( DofType v, DofType w ) { return (v - w < 1e-10); } ) )
          DUNE_THROW( Exception, "Local interpolation and basis function set are inconsistent" );
      }

      std::cout << std::endl;
      std::cout << "Phase II: Verifying that each DoF of the discrete function containts its global coordinate..." << std::endl;
      ConstLocalFunction< DiscreteFunctionType > uLocal( u );
      for( const auto &entity : space )
      {
        const auto &interpolation = space.interpolation( entity );

        idLocal.init( entity );
        vLocal.init( entity );
        interpolation( idLocal, vLocal.localDofVector() );

        uLocal.init( entity );
        wLocal.init( entity );
        interpolation( uLocal, wLocal.localDofVector() );

        if( !std::equal( vLocal.localDofVector().begin(), vLocal.localDofVector().end(), wLocal.localDofVector().begin(), [] ( DofType v, DofType w ) { return (v - w < 1e-10); } ) )
          DUNE_THROW( Exception, "Inconsistent Mapper" );
      }
    }

  } // namespace Fem

} // namespace Dune
