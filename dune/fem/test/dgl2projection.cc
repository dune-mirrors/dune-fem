#undef NDEBUG

#include <config.h>
#include <iostream>

#if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS && defined USE_IDGRIDPART
// do nothing in this test if experimental_grid_extension is not activated and it is tested for IDGridPart
int main () { return 0;}
#else

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/lpnorm.hh>

#if defined USE_BLOCKVECTORFUNCTION
#include <dune/fem/function/blockvectorfunction.hh>
#elif defined USE_VECTORFUNCTION
#include <dune/common/dynvector.hh>
#include <dune/fem/function/vectorfunction.hh>
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#elif HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
#include <dune/fem/function/petscdiscretefunction.hh>
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
  #if defined USE_DOFTYPE_INT && not defined USE_COMPLEX   // note: no conversion from int to std::complex implemented in std
        Dune :: DynamicVector< int >
  #else
        Dune :: DynamicVector< FunctionSpaceType :: RangeFieldType >
  #endif
  > >  DiscreteFunctionType;
#elif defined USE_BLOCKVECTORDISCRETEFUNCTION
typedef Dune::Fem::ReferenceBlockVector< FunctionSpaceType::RangeFieldType, DiscreteFunctionSpaceType::localBlockSize >
  BlockVectorType;
typedef Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType >
  DiscreteFunctionType;
#elif HAVE_PETSC && defined USE_PETSCDISCRETEFUNCTION
typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType >
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

typedef Fem::DiscreteFunctionTraits< DiscreteFunctionType > Traits;
typedef typename Traits::DofType DofType;

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
    interpolate( exactSolution, solution );

    LPNorm< GridPartType > lpnorm( gridPart, 2.0 );
    L2Norm< GridPartType > l2norm( gridPart );
    L1Norm< GridPartType > l1norm( gridPart );
    H1Norm< GridPartType > h1norm( gridPart );

    // weighted norm stuff
    WeightFunctionType weightFunctionExact;
    typedef GridFunctionAdapter< WeightFunctionType, GridPartType>
      DiscreteWeightFunctionType;

    typedef DiscreteFunctionType::RangeFieldType RangeFieldType ;

    DiscreteWeightFunctionType weightFunction( "weight", weightFunctionExact, gridPart );

    WeightedLPNorm< DiscreteWeightFunctionType > wLpnorm( weightFunction, 2.0 );
    WeightedL2Norm< DiscreteWeightFunctionType > wL2norm( weightFunction );

    // check all norms
    {
      // check lp norm
      RangeFieldType lperror  = lpnorm.distance( exactSolution, solution );
      RangeFieldType lperror2 = lpnorm.distance( solution, exactSolution );
      assert( std::abs( lperror - lperror2 ) < 1e-10 );

      // check l2 norm
      RangeFieldType error  = l2norm.distance( exactSolution, solution );
      RangeFieldType error2 = l2norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );

      // compare lp(p=2) and l2 norm
      assert( std::abs( lperror - error ) < 1e-10 );
    }

    // check l1 norm
    {
      RangeFieldType error  = l1norm.distance( exactSolution, solution );
      RangeFieldType error2 = l1norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );
    }

    // check h1 norm
    {
      RangeFieldType error  = h1norm.distance( exactSolution, solution );
      RangeFieldType error2 = h1norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );
    }

    // check weighted lp norm
    {
      // check lp norm
      RangeFieldType lperror  = wLpnorm.distance( exactSolution, solution );
      RangeFieldType lperror2 = wLpnorm.distance( solution, exactSolution );
      assert( std::abs( lperror - lperror2 ) < 1e-10 );

      // check l2 norm
      RangeFieldType error  = wL2norm.distance( exactSolution, solution );
      RangeFieldType error2 = wL2norm.distance( solution, exactSolution );
      assert( std::abs( error - error2 ) < 1e-10 );

      // compare lp(p=2) and l2 norm
      assert( std::abs( lperror - error ) < 1e-10 );
    }

    return 0;
  }
  catch( const Exception &e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
#endif // #if not DUNE_GRID_EXPERIMENTAL_GRID_EXTENSIONS && defined USE_IDGRIDPART
