#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/hierarchicgridpart.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#if defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#elif defined USE_VECTORFUNCTION
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction.hh>
#elif defined USE_ATTACHEDFUNCTION
#include <dune/fem/function/attachedfunction/function.hh>
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#elif defined USE_COMBINEDFUNCTION
#undef DIMRANGE 
#define DIMRANGE 1
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/combinedfunction.hh>
#else
#include <dune/fem/function/adaptivefunction.hh>
#endif

#if defined  USE_FILTEREDGRID 
#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#endif
#if defined  USE_IDGRIDPART 
#include <dune/fem/gridpart/idgridpart.hh>
#endif

#include "testgrid.hh"
#include "dfspace.hh"
#include "exactsolution.hh"
#include "weightfunction.hh"

using namespace Dune;
using namespace Fem;

typedef GridSelector::GridType MyGridType;
// typedef HierarchicGridPart< MyGridType >  ContainedGridPartType;
typedef Fem :: DGAdaptiveLeafGridPart< MyGridType > ContainedGridPartType;
//typedef IdBasedLeafGridPart< MyGridType > ContainedGridPartType;

// use filtered grid for testing 
#if defined  USE_FILTEREDGRID 
  typedef Fem :: RadialFilter< double, MyGridType :: dimensionworld > FilterImplType;
  typedef Fem :: BasicFilterWrapper< ContainedGridPartType, FilterImplType > FilterType ;
  typedef Fem :: FilteredGridPart<ContainedGridPartType, FilterType, true > GridPartType;
#elif defined  USE_IDGRIDPART
  typedef Fem:: IdGridPart< ContainedGridPartType > GridPartType;
#else 
  typedef ContainedGridPartType  GridPartType ;
#endif

typedef TestFunctionSpace FunctionSpaceType;
typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;

#if defined USE_BLOCKVECTORFUNCTION
typedef Fem :: ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_VECTORFUNCTION
typedef Fem :: ManagedDiscreteFunction
  < Fem :: VectorDiscreteFunction
    < DiscreteFunctionSpaceType,
  #if defined USE_DOFTYPE_INT
        Fem :: DynamicVector< int > 
  #else
        Fem :: DynamicVector< FunctionSpaceType :: RangeFieldType > 
  #endif
  > >  DiscreteFunctionType;
#elif defined USE_ATTACHEDFUNCTION
typedef Fem :: AttachedDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
typedef Dune::Fem::ReferenceBlockVector< double, DiscreteFunctionSpaceType::localBlockSize > 
  BlockVectorType;
typedef Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType > 
  DiscreteFunctionType;
#elif defined USE_COMBINEDFUNCTION
typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  ContainedDiscreteFunctionType;
typedef Fem :: CombinedDiscreteFunction< ContainedDiscreteFunctionType, DIMRANGE >
  DiscreteFunctionType;
#else
typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
  DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;
typedef Fem :: FunctionSpace< double, double, GridSelector::dimworld, 1 > WeightFunctionSpaceType;
typedef WeightFunction< WeightFunctionSpaceType > WeightFunctionType;

// main program 
int main(int argc, char ** argv) 
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );

#ifdef  USE_FILTEREDGRID 
    ContainedGridPartType containedGridPart( grid );
    FilterType filter( containedGridPart );
    GridPartType gridPart( containedGridPart, filter );
#else 
    GridPartType gridPart( grid );
#endif

    // add check for grid width 
    std::cout << "Grid width: " 
      << GridWidth :: calcGridWidth( gridPart ) << std::endl; 

    DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
    ExactSolutionType exactSolution;
    DiscreteFunctionType solution( "solution", discreteFunctionSpace );
    solution.clear();
  
    // perform the L2Projection
    Fem :: L2Projection< ExactSolutionType, DiscreteFunctionType > dgl2; 
    dgl2( exactSolution, solution );

    LPNorm< GridPartType > lpnorm( gridPart, 2.0 );
    L2Norm< GridPartType > l2norm( gridPart );
    L1Norm< GridPartType > l1norm( gridPart );
    H1Norm< GridPartType > h1norm( gridPart );

    // weighted norm stuff
    WeightFunctionType weightFunctionExact;
    typedef GridFunctionAdapter< WeightFunctionType, GridPartType>
      DiscreteWeightFunctionType;
      
    DiscreteWeightFunctionType weightFunction( "weight", weightFunctionExact, gridPart );

    WeightedLPNorm< DiscreteWeightFunctionType > wLpnorm( weightFunction, 2.0 );
    WeightedL2Norm< DiscreteWeightFunctionType > wL2norm( weightFunction );

    // check all norms
    {
      // check lp norm 
      double lperror  = lpnorm.distance( exactSolution, solution );
      double lperror2 = lpnorm.distance( solution, exactSolution );
      assert( std::abs( lperror - lperror2 ) < 1e-10 );

      // check l2 norm 
      double error  = l2norm.distance( exactSolution, solution );
      double error2 = l2norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );

      // compare lp(p=2) and l2 norm 
      assert( std::abs( lperror - error ) < 1e-10 );
    }

    // check l1 norm 
    {
      double error  = l1norm.distance( exactSolution, solution );
      double error2 = l1norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );
    }

    // check h1 norm 
    {
      double error  = h1norm.distance( exactSolution, solution );
      double error2 = h1norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );
    }

    // check weighted lp norm
    {
      // check lp norm 
      double lperror  = wLpnorm.distance( exactSolution, solution );
      double lperror2 = wLpnorm.distance( solution, exactSolution );
      assert( std::abs( lperror - lperror2 ) < 1e-10 );

      // check l2 norm 
      double error  = wL2norm.distance( exactSolution, solution );
      double error2 = wL2norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );

      // compare lp(p=2) and l2 norm 
      assert( std::abs( lperror - error ) < 1e-10 );
    }

    return 0;
  }
  catch( Exception e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
