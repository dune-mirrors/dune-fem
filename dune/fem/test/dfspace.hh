#ifndef DUNE_FEM_TEST_DFSPACE_HH
#define DUNE_FEM_TEST_DFSPACE_HH

#if !HAVE_GRIDTYPE
#define DEFAULT_GRID
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>
#endif

#ifndef GRIDDIM
#define GRIDDIM 2
#endif

#include <complex>

#include <dune/fem/misc/double.hh>
#include <dune/fem/space/common/functionspace.hh>

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
#define POLORDER 2
#endif

#ifndef DIMRANGE
#define DIMRANGE 4
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

#ifdef DEFAULT_GRID
#if not defined USE_COMPLEX
  typedef FunctionSpace< double, RangeFieldType, GRIDDIM, DIMRANGE >
    TestFunctionSpace;
#else
  typedef FunctionSpace< double, std::complex<RangeFieldType>, GRIDDIM, DIMRANGE >
    TestFunctionSpace;
#endif
#else
#if not defined USE_COMPLEX
  typedef FunctionSpace< double, RangeFieldType, GridSelector::dimworld, DIMRANGE >
    TestFunctionSpace;
#else
  typedef FunctionSpace< double, std::complex<RangeFieldType>, GridSelector::dimworld, DIMRANGE >
    TestFunctionSpace;
#endif
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
}

#endif
