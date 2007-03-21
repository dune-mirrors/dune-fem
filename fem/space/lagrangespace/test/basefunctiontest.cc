#include "basefunctiontest.hh"

#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/quadrature/cachequad.hh>

#include <dune/fem/space/lagrangespace/lagrangepoints.hh>

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
  
    // check polynomial order 2
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 >
      TwoSpaceType;
    {
      std :: cout << "Quadratic Base Functions" << std :: endl;
      TwoSpaceType space( gridPart );
      checkLagrangeBase( space );
    }
  }
  
  template< class SpaceType > 
  void LagrangeBase_Test :: checkLagrangeBase( const SpaceType &space )
  {
    typedef typename SpaceType :: IteratorType IteratorType; 
    typedef typename SpaceType :: BaseFunctionSetType
      BaseFunctionSetType;
    typedef typename SpaceType :: GridPartType GridPartType;
    typedef typename SpaceType :: DomainType DomainType; 
    typedef typename SpaceType :: RangeType RangeType; 
    //typedef typename SpaceType :: GridType :: template Codim< 0 > :: Geometry
    //  GeometryType;

    int errors = 0;
    
    IteratorType end = space.end();
    for(IteratorType it = space.begin(); it != end; ++it)
    {
      const BaseFunctionSetType& baseSet = space.baseFunctionSet( *it );
      const int numBaseFct = baseSet.numBaseFunctions();

      //const GeometryType& geo = it->geometry();
     
      LagrangeQuadrature< GridPartType, 0 >
        lagrangePoints( *it, space.order() );

      for( int i = 0; i < numBaseFct; ++i ) 
      {
        RangeType phi( 0.0 );
       
        const DomainType& x = lagrangePoints.point( i );
        
        // eval on lagrange point 
        baseSet.eval( i , x , phi ); 
        if( std :: abs( phi[ 0 ] - 1.0 ) >= 1e-10 ) {
          std :: cout << "Base function " << i << " failed at " << x << " (" << phi[ 0 ] << " != 1)!" << std :: endl;
          errors++;
        }
        
        for( int j = 0; j < numBaseFct; ++j ) 
        {
          if( i == j )
            continue;

          // eval on lagrange point 
          baseSet.eval( j , x, phi ); 
          if( std :: abs( phi[ 0 ] ) >= 1e-10 ) {
            std :: cout << "Base function " << j << " failed at " << x << " (" << phi[ 0 ] << " != 0)!" << std :: endl;
            errors++;
          }
        }
      }
    }

    assert( errors == 0 );
  }

} // end namespace Dune
