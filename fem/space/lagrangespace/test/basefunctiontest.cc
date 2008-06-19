#include "basefunctiontest.hh"

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>

namespace Dune {

  void LagrangeBase_Test::run() 
  {
    testBaseFunctions();
  }

  void LagrangeBase_Test::testBaseFunctions() 
  {
    GridPtr< GridType > gridPtr( gridFile_ );
    GridType& grid = *gridPtr;

    typedef LeafGridPart< GridType > GridPartType; 
    GridPartType gridPart( grid );

    typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;
    
    // check polynomial order 1 
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >
      OneSpaceType;
    {
      std :: cout << "Linear Base Functions" << std :: endl;
      OneSpaceType space( gridPart);
      checkLagrangeBase( space );
    }
 
    #ifdef TEST_SECOND_ORDER
    // check polynomial order 2
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 >
      TwoSpaceType;
    {
      std :: cout << "Quadratic Base Functions" << std :: endl;
      TwoSpaceType space( gridPart );
      checkLagrangeBase( space );
    }
    #endif

    #ifdef TEST_THIRD_ORDER
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 3 >
      ThreeSpaceType;
    {
      std :: cout << "Cubic Base Functions" << std :: endl;
      ThreeSpaceType space( gridPart );
      checkLagrangeBase( space );
    }
    #endif

  }
  
  template< class SpaceType > 
  void LagrangeBase_Test :: checkLagrangeBase( const SpaceType &space )
  {
    typedef typename SpaceType :: IteratorType IteratorType; 
    typedef typename SpaceType :: BaseFunctionSetType  BaseFunctionSetType;
    typedef typename SpaceType :: LagrangePointSetType LagrangePointSetType;
    typedef typename SpaceType :: GridPartType GridPartType;
    typedef typename SpaceType :: DomainType DomainType; 
    typedef typename SpaceType :: RangeType RangeType; 

    int errors = 0;
    
    IteratorType end = space.end();
    for(IteratorType it = space.begin(); it != end; ++it)
    {
      const BaseFunctionSetType& baseSet = space.baseFunctionSet( *it );
      const int numBaseFct = baseSet.numBaseFunctions();

      const LagrangePointSetType& pointSet = space.lagrangePointSet( *it );
      const int numPoints = pointSet.size();
     
      for( int i = 0; i < numPoints; ++i ) 
      {
        RangeType phi( 0.0 );
       
        const DomainType& x = pointSet.point( i );
       
        // evaluate on lagrange point 
        baseSet.evaluate( i , x , phi ); 
        if( std :: abs( phi[ 0 ] - 1.0 ) >= 1e-10 )
        {
          std :: cout << "Base function " << i << " failed at " << x
                      << " (" << phi[ 0 ] << " != 1)!" << std :: endl;
          ++errors;
        }
        
        for( int j = 0; j < numBaseFct; ++j ) 
        {
          if( i == j )
            continue;

          // evaluate on lagrange point 
          baseSet.evaluate( j , x, phi ); 
          if( std :: abs( phi[ 0 ] ) >= 1e-10 )
          {
            std :: cout << "Base function " << j << " failed at " << x
                        << " (" << phi[ 0 ] << " != 0)!" << std :: endl;
            ++errors;
          }
        }
      }
    }

    assert( errors == 0 );
  }

} // end namespace Dune
