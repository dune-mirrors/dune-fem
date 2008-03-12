#ifndef DUNE_FEM_TEST_DFSPACE_HH
#define DUNE_FEM_TEST_DFSPACE_HH

#include <dune/fem/misc/double.hh>

#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/common/dofstorage.hh>
#if defined USE_COMBINEDSPACE
#  include <dune/fem/space/combinedspace.hh>
#endif
#if defined USE_FVSPACE
#  include <dune/fem/space/fvspace.hh>
#endif


#ifndef POLORDER
#  define POLORDER 1
#endif

#ifndef DIMRANGE
#  define DIMRANGE 2
#endif


namespace Dune
{

  typedef FunctionSpace< double, Double, dimworld, DIMRANGE >
    TestFunctionSpace;


  template< class GridPart >
  struct TestDiscreteFunctionSpaceTraits
  {
    enum { polOrder = POLORDER };

    typedef GridPart GridPartType;

    typedef TestFunctionSpace FunctionSpaceType;
    
    typedef DiscontinuousGalerkinSpace
      < FunctionSpaceType :: ScalarFunctionSpaceType, GridPartType, polOrder >
      SingleDiscreteFunctionSpaceType;
    enum { dimRange = FunctionSpaceType :: dimRange };

#ifdef USE_VARIABLEBASE
    static const DofStoragePolicy policy = VariableBased;
#else
    static const DofStoragePolicy policy = PointBased;
#endif

#if defined USE_COMBINEDSPACE
    typedef CombinedSpace< SingleDiscreteFunctionSpaceType, dimRange, policy >
      DiscreteFunctionSpaceType;
#elif defined USE_FVSPACE
    typedef FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0 >
      DiscreteFunctionSpaceType;
#else
    typedef DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;
#endif
  };


  template< class GridPart >
  class TestDiscreteFunctionSpace
  : public TestDiscreteFunctionSpaceTraits< GridPart > :: DiscreteFunctionSpaceType
  {
  public:
    typedef GridPart GridPartType;

  private:
    typedef typename TestDiscreteFunctionSpaceTraits< GridPartType >
      :: DiscreteFunctionSpaceType
      BaseType;

  public:
    TestDiscreteFunctionSpace ( GridPartType &gridPart )
    : BaseType( gridPart )
    {}
  };
    
}

#endif
