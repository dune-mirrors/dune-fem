#include "mappertest.hh"

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>

namespace Dune {

  void LagrangeMapper_Test :: run()
  {
    GridPtr< GridType > gridPtr( gridFile_ );
    GridType& grid = *gridPtr;
    grid.globalRefine( 2 );
    GridPartType gridPart( grid );

    typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

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
  }


  
  template< class SpaceType >
  void LagrangeMapper_Test :: checkDiscreteFunction( const SpaceType &space )
  {
    typedef AdaptiveDiscreteFunction< SpaceType > DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType :: LocalFunctionType LocalFunctionType;
    typedef typename DiscreteFunctionType :: DofIteratorType DofIteratorType;
    typedef typename SpaceType :: IteratorType IteratorType;

    DiscreteFunctionType u( "u", space );
    u.clear();

    int errors = 0;

    IteratorType eit = space.end();
    for( IteratorType it = space.begin(); it != eit; ++it )
    {
      LocalFunctionType ulocal = u.localFunction( *it );
      
      const int numDofs = ulocal.numDofs();
      std :: cout << numDofs << "  ";
      for( int i = 0; i < numDofs; ++i )
        ulocal[ i ] += 1.0;
    }

    DofIteratorType deit = u.dend();
    for( DofIteratorType dit = u.dbegin(); dit != deit; ++dit ) {
      if( *dit < 0.5 )
        errors++;
      std :: cout << *dit << "  ";
    }
    std :: cout << std :: endl;

    assert( errors == 0 );
  }

    
} // End namespace Dune
