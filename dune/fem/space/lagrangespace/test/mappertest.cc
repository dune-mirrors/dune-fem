#include "mappertest.hh"

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune
{

  template< class Grid >
  void LagrangeMapper_Test< Grid > :: run()
  {
    GridPtr< GridType > gridPtr( gridFile_ );
    GridType& grid = *gridPtr;
    //grid.globalRefine( 2 );
    GridPartType gridPart( grid );

    typedef FunctionSpace< double, double, dimworld, dimworld > FunctionSpaceType;

    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >
      OneSpaceType;
    {
      std :: cout << "Testing linear function." << std :: endl;
      OneSpaceType space( gridPart );
      checkDiscreteFunction( space );
    }

    #ifdef TEST_SECOND_ORDER
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 >
      TwoSpaceType;
    {
      std :: cout << "Testing quadratic function." << std :: endl;
      TwoSpaceType space( gridPart );
      checkDiscreteFunction( space );
    }
    #endif

    #ifdef TEST_THIRD_ORDER
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 3 >
      ThreeSpaceType;
    {
      std :: cout << "Testing cubic function." << std :: endl;
      ThreeSpaceType space( gridPart );
      checkDiscreteFunction( space );
    }
    #endif
  }



  template< class Grid >
  template< class SpaceType >
  void LagrangeMapper_Test< Grid >
    :: checkDiscreteFunction( const SpaceType &space )
  {
    typedef AdaptiveDiscreteFunction< SpaceType > DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
    typedef typename SpaceType :: LagrangePointSetType LagrangePointSetType;
    typedef typename SpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;

#if 0
    typedef typename SpaceType :: GridType GridType;
    for( unsigned int i = 0; i <= GridType :: dimension; ++i )
      std :: cout << "size of codimension " << i << ": "
                  << space.grid().size( i ) << std :: endl;
#endif

    std :: cout << "size of space: " << space.size() << " (= "
                << (space.size() / dimworld) << " * " << dimworld << ")"
                << std :: endl;

    DiscreteFunctionType u( "u", space );
    u.clear();

    int errors = 0;

    std :: cout << std :: endl << "Phase I: "
                << "Setting each DoF of a discrete function to its global "
                << "coordinate..." << std :: endl;

    const IteratorType eit = space.end();
    for( IteratorType it = space.begin(); it != eit; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      const LagrangePointSetType &lagrangePoints
        = space.lagrangePointSet( entity );
      const int numLPoints = lagrangePoints.nop();

      LocalFunctionType ulocal = u.localFunction( entity );

      assert( numLPoints * dimworld == ulocal.numDofs() );
      for( int i = 0; i < numLPoints; ++i )
      {
        const FieldVector< double, dimworld > &lpoint
          = lagrangePoints.point( i );
        FieldVector< double, dimworld > x
          = geometry.global( lpoint );

        for( int j = 0; j < dimworld; ++j )
          ulocal[ i * dimworld + j ] = x[ j ];

        FieldVector< double, dimworld > y( 0 );
        ulocal.evaluate( lagrangePoints[ i ], y );
        if( (y - x).two_norm() > 1e-10 )
        {
          std :: cout << "point " << i << " ( " << lpoint << " ): "
                      << x << " != " << y << std :: endl;
          ++errors;
        }
      }
    }

    std :: cout << std :: endl << "Phase II: "
                << "Verifying that each DoF of the discrete function "
                << "containts its global" << std :: endl
                << "          coordinate..." << std :: endl;
    for( IteratorType it = space.begin(); it != eit; ++it )
    {
      const EntityType &entity = *it;
      const GeometryType &geometry = entity.geometry();

      const LagrangePointSetType &lagrangePoints
        = space.lagrangePointSet( entity );
      const int numLPoints = lagrangePoints.nop();

      LocalFunctionType ulocal = u.localFunction( entity );
      
      assert( numLPoints * dimworld == ulocal.numDofs() );
      for( int i = 0; i < numLPoints; ++i )
      {
        const FieldVector< double, dimworld > &lpoint
          = lagrangePoints.point( i );
        FieldVector< double, dimworld > x
          = geometry.global( lpoint );

        FieldVector< double, dimworld > y( 0 );
        ulocal.evaluate( lagrangePoints[ i ], y );
        if( (y - x).two_norm() > 1e-10 )
        {
          std :: cout << "point " << i << " ( " << lpoint << " ): "
                      << x << " != " << y << std :: endl;
          ++errors;
        }
      }
    }

    assert( errors == 0 );
  }
    
} // End namespace Dune
