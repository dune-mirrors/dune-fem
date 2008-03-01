#include "mappertest.hh"

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune
{

  void LagrangeMapper_Test :: run()
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


  
  template< class SpaceType >
  void LagrangeMapper_Test :: checkDiscreteFunction( const SpaceType &space )
  {
    typedef AdaptiveDiscreteFunction< SpaceType > DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
    typedef typename SpaceType :: LagrangePointSetType LagrangePointSetType;
    typedef typename SpaceType :: IteratorType IteratorType;
    typedef typename IteratorType :: Entity EntityType;
    typedef typename EntityType :: Geometry GeometryType;

    typedef typename SpaceType :: GridType GridType;
    for( unsigned int i = 0; i <= GridType :: dimension; ++i )
      std :: cout << "size of codimension " << i << ": "
                  << space.grid().size( i ) << std :: endl;

    std :: cout << "size of space: " << space.size() << std :: endl;

    DiscreteFunctionType u( "u", space );
    u.clear();

    int errors = 0;

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
        FieldVector< double, dimworld > x
          = geometry.global( lagrangePoints.point( i ) );
        for( int j = 0; j < dimworld; ++j )
          ulocal[ i * dimworld + j ] = x[ j ];
      }
    }

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
        FieldVector< double, dimworld > x
          = geometry.global( lagrangePoints.point( i ) );
        FieldVector< double, dimworld > y;
        ulocal.evaluate( lagrangePoints[ i ], y );
        if( (y - x).two_norm() > 1e-6 )
        {
          std :: cout << x << " != " << y << std :: endl;
          ++errors;
        }
      }
    }

    assert( errors == 0 );
  }

    
} // End namespace Dune
