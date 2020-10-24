// undef NDEBUG so we can use assert every time.
#undef NDEBUG

#include <config.h>
#include <iostream>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fem/misc/l1norm.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/misc/lpnorm.hh>

#include <dune/common/dynvector.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/combinedfunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/petscdiscretefunction.hh>

#include <dune/fem/gridpart/filter/radialfilter.hh>
#include <dune/fem/gridpart/filter/basicfilterwrapper.hh>
#include <dune/fem/gridpart/filteredgridpart.hh>
#include <dune/fem/gridpart/idgridpart.hh>


#include <complex>
#include <dune/fem/misc/double.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/common/dofstorage.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/finitevolume.hh>

#include "testgrid.hh"
#include "dfspace.hh"
#include "exactsolution.hh"
#include "weightfunction.hh"

using namespace Dune;
using namespace Fem;

typedef GridSelector::GridType MyGridType;
typedef Fem :: DGAdaptiveLeafGridPart< MyGridType > ContainedGridPartType;

template <class DiscreteFunctionType>
void algorithm( DiscreteFunctionType& solution )
{
  typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType :: GridPartType GridPartType;

  // just to check that it is there
  typedef Fem::DiscreteFunctionTraits< DiscreteFunctionType > Traits;
  typedef typename Traits::DofType DofType;
  DofType val = 0;
  val = val + val;

  typedef ExactSolution< typename DiscreteFunctionSpaceType::FunctionSpaceType > ExactSolutionType;
  typedef Fem :: FunctionSpace< double, double, GridPartType::dimensionworld, 1 > WeightFunctionSpaceType;
  typedef WeightFunction< WeightFunctionSpaceType > WeightFunctionType;

  const GridPartType& gridPart = solution.space().gridPart();

  ExactSolutionType exactSolution;
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

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType ;

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
}

// default run algorithm
template <class GridPartType, class DiscreteFunctionSpaceType, class DiscreteFunctionType>
void run(MyGridType& grid)
{
  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  algorithm( solution );
}

// main program
int main(int argc, char ** argv)
{
  MPIManager :: initialize( argc, argv );
  try
  {
    MyGridType &grid = TestGrid :: grid();
    const int step = TestGrid :: refineStepsForHalf();
    grid.globalRefine( 2*step );

    // AdaptiveDiscreteFunction
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    // AdaptiveDiscreteFunction + FV space
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef FiniteVolumeSpace< TestFunctionSpace, GridPartType, 0 > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    // AdaptiveDiscreteFunction + DG HierarchicLegendre
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef HierarchicLegendreDiscontinuousGalerkinSpace<
        TestFunctionSpace, GridPartType, 2 > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    // ISTLBlockVectorDiscreteFunction
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    // ManagedDiscreteFunction + int
    {
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Dune :: DynamicVector< int >  VectorType;
      typedef Fem :: ManagedDiscreteFunction< Fem :: VectorDiscreteFunction < DiscreteFunctionSpaceType, VectorType> > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    // ManagedDiscreteFunction + RangeFieldType
    {
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Dune :: DynamicVector< typename DiscreteFunctionSpaceType :: RangeFieldType >  VectorType;
      typedef Fem :: ManagedDiscreteFunction< Fem :: VectorDiscreteFunction < DiscreteFunctionSpaceType, VectorType> > DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

    /*
    // CombinedDiscreteFunction
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType >
        ContainedDiscreteFunctionType;
      typedef Fem :: CombinedDiscreteFunction< ContainedDiscreteFunctionType, 5 >
        DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }
    */

    // BlockVectorDiscreteFunction
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Dune::Fem::ReferenceBlockVector< typename DiscreteFunctionSpaceType::RangeFieldType, DiscreteFunctionSpaceType::localBlockSize >
        BlockVectorType;
      typedef Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType >
        DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }

#if HAVE_PETSC
    // PetscDiscreteFunction
    {
      // use standard grid part
      typedef ContainedGridPartType  GridPartType ;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType >
        DiscreteFunctionType;
      run< GridPartType, DiscreteFunctionSpaceType, DiscreteFunctionType > ( grid );
    }
#endif

    // AdaptiveDiscreteFunction + FilteredGridPart
    {
      // use standard grid part
      typedef Fem :: RadialFilter< double, MyGridType :: dimensionworld > FilterImplType;
      typedef Fem :: BasicFilterWrapper< ContainedGridPartType, FilterImplType > FilterType ;
      typedef Fem :: FilteredGridPart<ContainedGridPartType, FilterType, true > GridPartType;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

      ContainedGridPartType containedGridPart( grid );
      FilterType filter( containedGridPart );
      GridPartType gridPart( containedGridPart, filter );

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      DiscreteFunctionType solution( "solution", discreteFunctionSpace );
      algorithm( solution );
    }

    // ISTLBlockVectorDiscreteFunction + FilteredGridPart
    {
      // use standard grid part
      typedef Fem :: RadialFilter< double, MyGridType :: dimensionworld > FilterImplType;
      typedef Fem :: BasicFilterWrapper< ContainedGridPartType, FilterImplType > FilterType ;
      typedef Fem :: FilteredGridPart<ContainedGridPartType, FilterType, true > GridPartType;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

      ContainedGridPartType containedGridPart( grid );
      FilterType filter( containedGridPart );
      GridPartType gridPart( containedGridPart, filter );

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      DiscreteFunctionType solution( "solution", discreteFunctionSpace );
      algorithm( solution );
    }

    // AdaptiveDiscreteFunction + IdGridPart
    {
      // use standard grid part
      typedef Fem:: IdGridPart< ContainedGridPartType > GridPartType;
      typedef TestDiscreteFunctionSpace< GridPartType > DiscreteFunctionSpaceType;
      typedef Fem :: AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

      ContainedGridPartType containedGridPart( grid );
      GridPartType gridPart( containedGridPart );

      DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
      DiscreteFunctionType solution( "solution", discreteFunctionSpace );
      algorithm( solution );
    }

    return 0;
  }
  catch( const Exception &e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
