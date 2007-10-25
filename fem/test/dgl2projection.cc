#include <config.h>
#include <iostream>

#include <dune/grid/common/gridpart.hh>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/dgspace.hh>
#ifdef USE_COMBINEDSPACE
#include <dune/fem/space/combinedspace.hh>
#endif
#ifdef USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#endif

#include "testgrid.hh"
#include "dgl2projection.hh"
#include "exactsolution.hh"

using namespace Dune;

// polynom approximation order of quadratures, 
// at least poolynom order of basis functions 
#ifdef POLORDER
  const int polOrder = POLORDER;
#else
  const int polOrder = 1;
#endif

typedef HierarchicGridPart< GridType > GridPartType;

#ifdef DIMRANGE
  const int dimRange = DIMRANGE;
#else
  const int dimRange = 2;
#endif

typedef FunctionSpace< double, Double, dimworld, dimRange > FunctionSpaceType;
#ifdef USE_COMBINEDSPACE
//typedef FunctionSpace < double , Double, dimworld, 1 > SingleFunctionSpace;
typedef DiscontinuousGalerkinSpace
  < FunctionSpaceType :: ScalarFunctionSpaceType, GridPartType, polOrder >
  SingleDiscreteFunctionSpaceType;
#ifdef USE_VARIABLEBASE
typedef CombinedSpace< SingleDiscreteFunctionSpaceType, dimRange, VariableBased >
   DiscreteFunctionSpaceType;
#else
typedef CombinedSpace< SingleDiscreteFunctionSpaceType, dimRange, PointBased >
   DiscreteFunctionSpaceType;
#endif
#else
typedef DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;
#endif
  
#ifdef USE_BLOCKVECTORFUNCTION
typedef BlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#else
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

int main ()
{
  try
  {
    GridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();

    grid.globalRefine( 2*step );

    GridPartType gridPart( grid );
    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType exactSolution( discreteFunctionSpace );
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();
  
    // perform the L2Projection
    DGL2Projection< DiscreteFunctionType > :: project( exactSolution, solution );

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
