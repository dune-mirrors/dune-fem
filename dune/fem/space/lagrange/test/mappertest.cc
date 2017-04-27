#include "mappertest.hh"

#include <dune/fem/function/adaptivefunction.hh>
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

      typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
      typedef typename SpaceType::LagrangePointSetType LagrangePointSetType;
      typedef typename SpaceType::IteratorType::Entity::Geometry GeometryType;
      static const int dimworld = SpaceType::GridPartType::GridType::dimensionworld;


      std::cout << "size of space: " << space.size() << " (= " << (space.size() / dimworld) << " * " << dimworld << ")" << std::endl;

      DiscreteFunctionType u( "u", space );
      u.clear();

      int errors = 0;

      std::cout << std::endl;
      std::cout << "Phase I: Setting each DoF of a discrete function to its global coordinate..." << std::endl;

      for( const auto &entity : space )
      {
        const GeometryType &geometry = entity.geometry();

        const LagrangePointSetType &lagrangePoints = space.lagrangePointSet( entity );
        const int numLPoints = lagrangePoints.nop();

        LocalFunctionType ulocal = u.localFunction( entity );

        assert( numLPoints * dimworld == ulocal.numDofs() );
        for( int i = 0; i < numLPoints; ++i )
        {
          const FieldVector< double, dimworld > &lpoint = lagrangePoints.point( i );
          FieldVector< double, dimworld > x = geometry.global( lpoint );

          for( int j = 0; j < dimworld; ++j )
            ulocal[ i * dimworld + j ] = x[ j ];

          FieldVector< double, dimworld > y( 0 );
          ulocal.evaluate( lagrangePoints[ i ], y );
          if( (y - x).two_norm() > 1e-10 )
          {
            std::cout << "point " << i << " ( " << lpoint << " ): " << x << " != " << y << std::endl;
            ++errors;
          }
        }
      }

      std::cout << std::endl << "Phase II: Verifying that each DoF of the discrete function containts its global coordinate..." << std::endl;
      for( const auto &entity : space )
      {
        const GeometryType &geometry = entity.geometry();

        const LagrangePointSetType &lagrangePoints = space.lagrangePointSet( entity );
        const int numLPoints = lagrangePoints.nop();

        LocalFunctionType ulocal = u.localFunction( entity );

        assert( numLPoints * dimworld == ulocal.numDofs() );
        for( int i = 0; i < numLPoints; ++i )
        {
          const FieldVector< double, dimworld > &lpoint = lagrangePoints.point( i );
          FieldVector< double, dimworld > x = geometry.global( lpoint );

          FieldVector< double, dimworld > y( 0 );
          ulocal.evaluate( lagrangePoints[ i ], y );
          if( (y - x).two_norm() > 1e-10 )
          {
            std::cout << "point " << i << " ( " << lpoint << " ): " << x << " != " << y << std::endl;
            ++errors;
          }
        }
      }

      assert( errors == 0 );
    }

  } // namespace Fem

} // namespace Dune
