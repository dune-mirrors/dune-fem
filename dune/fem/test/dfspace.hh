#ifndef DUNE_FEM_TEST_DFSPACE_HH
#define DUNE_FEM_TEST_DFSPACE_HH

#include <complex>

#include <dune/fem/misc/double.hh>

#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/common/dofstorage.hh>
#if defined USE_COMBINEDSPACE
#  include <dune/fem/space/combinedspace.hh>
#endif
#if defined USE_LAGRANGESPACE
#  include <dune/fem/space/lagrange.hh>
#endif
#if defined USE_FVSPACE
#  include <dune/fem/space/finitevolume.hh>
#endif


#ifndef POLORDER
#  define POLORDER 2
#  define DEFAULTPOLORDER 1
#endif

#ifndef DIMRANGE
#  define DIMRANGE 1
#  define DEFAULTDIMRANGE 1
#endif


namespace Dune
{

  namespace Fem
  {

// complex range type not yet working with petsc
#if defined USE_PETSCDISCRETEFUNCTION && defined USE_COMPLEX
#warning PETSC does not work yet with complex numbers - undefining USE_COMPLEX for now
#undef USE_COMPLEX
#endif

//#ifdef COUNT_FLOPS
//  typedef Dune::Fem::Double  RangeFieldType ;
//#else
  typedef double RangeFieldType;
//#endif

#if not defined USE_COMPLEX
  typedef FunctionSpace< double, RangeFieldType, GridSelector::dimworld, DIMRANGE >
    TestFunctionSpace;
#else
  typedef FunctionSpace< double, std::complex<RangeFieldType>, GridSelector::dimworld, DIMRANGE >
    TestFunctionSpace;
#endif

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
    typedef Fem::CombinedSpace< SingleDiscreteFunctionSpaceType, dimRange, policy >
      DiscreteFunctionSpaceType;
#elif defined USE_FVSPACE
    typedef FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0 >
      DiscreteFunctionSpaceType;
#elif defined USE_LAGRANGESPACE
#warning USING LAGRANGE
    typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;
#elif defined USE_LEGENDRESPACE
    typedef LegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;
#elif defined USE_HIERARCHICLEGENDRESPACE
    typedef HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;
#else
    typedef DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
      DiscreteFunctionSpaceType;
#endif
  };


  template< class GridPart >
  using TestDiscreteFunctionSpace = typename TestDiscreteFunctionSpaceTraits< GridPart > :: DiscreteFunctionSpaceType;

  }
}

#endif
