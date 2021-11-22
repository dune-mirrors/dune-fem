#include <config.h>

#define SHOW_INTERPOLATION 0
#define SHOW_RESTRICT_PROLONG 1

// to write out the data, set WRITE_DATA to 1
#define WRITE_DATA 0

// to use grape, set to WANT_GRAPE to 1
#ifndef WANT_GRAPE
#define WANT_GRAPE 0
#endif

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <sstream>
//#include <dune/common/stdstreams.cc>
#include <tuple>

#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#if HAVE_GRAPE
  #define USE_GRAPE WANT_GRAPE
#else
  #define USE_GRAPE 0
  #if WANT_GRAPE
    #warning "Grape was not found by configure."
  #endif
#endif

#if USE_GRAPE
  #include <dune/grid/io/visual/grapedatadisplay.hh>
#endif
#include <dune/fem/io/parameter.hh>

#include <dune/fem/test/testgrid.hh>

// Check for unhealthy grids
// -------------------------

// forward declaration of the real main method
int Main ( int argc, char **argv );

template< class Grid >
struct CheckGridEnabled
{
  typedef Grid GridType;

  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

  inline static int CallMain ( int argc, char **argv )
  {
    return Main( argc, argv );
  }
};

// disable YaspGrid
namespace Dune
{
  template< int dim, class CoordCont >
  class YaspGrid;
}

template< int dim >
struct CheckGridEnabled< Dune :: YaspGrid< dim > >
{
  typedef Dune :: YaspGrid< dim > GridType;

  typedef Dune :: Fem :: LeafGridPart< GridType > GridPartType;

  inline static int CallMain ( int argc, char **argv )
  {
    std :: cerr << "WARNING: Lagrange Adaptation test disabled, because YaspGrid sucks!"
                << std :: endl;
    return 0;
  }
};

int main ( int argc, char **argv )
{
  return CheckGridEnabled< Dune::GridSelector::GridType >::CallMain( argc, argv );
}

// use deprecated interface check unit ALUGrid returns entities instead of entity pointers
#ifndef DUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE
#define DUNE_GRID_CHECK_USE_DEPRECATED_ENTITY_AND_INTERSECTION_INTERFACE 1
#endif

#include <dune/grid/test/checkindexset.hh>
template <class GridPart>
void checkAdaptiveIndexSet( const GridPart& gridPart )
{
  // call check index set from the DUNE grid test suite
  Dune :: checkIndexSet( gridPart.grid(), gridPart, Dune :: dvverb );
}

using namespace Dune;
using namespace Fem;

// Function, we will interpolate
// -----------------------------

template< class FunctionSpace >
class ExactSolution
: public Fem::Function< FunctionSpace, ExactSolution< FunctionSpace > >
{
  typedef ExactSolution< FunctionSpace > ThisType;
  typedef Fem::Function< FunctionSpace, ThisType > BaseType;

public:
  typedef FunctionSpace FunctionSpaceType;

  typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType::DomainType DomainType;
  typedef typename FunctionSpaceType::RangeType RangeType;
  typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;

public:
  void evaluate ( const DomainType &x, RangeType &phi ) const
  {
    phi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      // phi[ 0 ] += x[ i ] * x[ i ];
      phi[ 0 ] *= sin( M_PI * x[ i ] );
  }

  void evaluate ( const DomainType &x, RangeFieldType t, RangeType &phi ) const
  {
    evaluate( x, phi );
  }

  void jacobian( const DomainType &x, JacobianRangeType &Dphi ) const
  {
    Dphi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      for( int j = 0; j < DomainType :: dimension; ++j )
        // Dphi[ 0 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
  }

  void jacobian( const DomainType &x, RangeFieldType t, JacobianRangeType &Dphi ) const
  {
    jacobian( x, Dphi );
  }
};



// Type Definitions
// ----------------

typedef Dune::GridSelector::GridType MyGridType;

typedef CheckGridEnabled< MyGridType >::GridPartType GridPartType;

//! type of the function space
typedef FunctionSpace< double, double, MyGridType::dimensionworld, 1 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder+1 >
  DiscreteFunctionSpaceTwoType;

//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceTwoType > DiscreteFunctionTwoType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef GridFunctionAdapter< ExactSolutionType, GridPartType >
  GridExactSolutionType;

//! type of the DoF manager
typedef DofManager< MyGridType > DofManagerType;

//! type of restrict-prolong operator
typedef RestrictProlongDefaultTuple<DiscreteFunctionType, DiscreteFunctionTwoType, DiscreteFunctionType>
RestrictProlongOperatorType;

//! type of the adaption manager
typedef AdaptationManager< MyGridType, RestrictProlongOperatorType >
  AdaptationManagerType;






template< class... DiscreteFunction >
void adapt( MyGridType &grid, std::tuple< DiscreteFunction &... > functionTuple, int step,
            const bool locallyAdaptive )
{
  DiscreteFunctionType& solution(std::get<0>(functionTuple));

  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();

  RestrictProlongOperatorType rp(functionTuple);
  RestrictProlongOperatorType rpTwo(std::get<0>(functionTuple),
                                    std::get<1>(functionTuple),
                                    std::get<2>(functionTuple));
  (void)rpTwo;

  AdaptationManagerType adaptationManager( grid, rp );

  std :: string message = (step < 0 ? "Coarsening..." : "Refining..." );
  const int mark = (step < 0 ? -1 : 1);
  const int count = std :: abs( step );

  for( int i = 0; i < count; ++i )
  {
    int numElements = grid.size( 0 );
    if( locallyAdaptive )
    {
      numElements /= 4;
      numElements = std::max( numElements, 1 );
    }

    int elementNumber = 0;
    for( const auto& entity : discreteFunctionSpace )
    {
      if( elementNumber < numElements )
        grid.mark( mark, entity );
      ++elementNumber;
    }

    // adapt grid
    adaptationManager.adapt();
  }
}



template< class... DiscreteFunction >
void algorithm ( GridPartType &gridPart,
                 std::tuple< DiscreteFunction &... > functionTuple,
                 int step,
                 int turn,
                 const bool locallyAdaptive )
{
  const unsigned int polOrder
    = DiscreteFunctionSpaceType :: maxPolynomialOrder + 1;

  DiscreteFunctionType &solution = std::get<0>(functionTuple);

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );

  interpolate( f, solution );
  double preL2error = l2norm.distance( f, solution );
  double preH1error = h1norm.distance( f, solution );

  std::cout << "Unknowns before adaptation: " << solution.space().size() << std::endl;
  std::cout << "L2 error before adaptation: " << preL2error << std::endl;
  std::cout << "H1 error before adaptation: " << preH1error << std::endl;

  DiscreteFunctionTwoType &second = std::get<1>(functionTuple);
  interpolate( f, second );

  adapt( gridPart.grid(), functionTuple, step, locallyAdaptive );

  // check index set for being consistent with the dune index set check
  checkAdaptiveIndexSet( gridPart );

  //double postL2error = l2norm.distance( fexact, solution );
  double postL2error = l2norm.distance( solution, fexact );
  double postH1error = h1norm.distance( f, solution );

  std::cout << "Unknowns after "
            << (step < 0 ? "restriction" : "prolongation")
            << ": " << solution.space().size() << std::endl;
  std::cout << "L2 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postL2error << std::endl;
  std::cout << "H1 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postH1error << std::endl;

  #if USE_GRAPE && SHOW_RESTRICT_PROLONG
    if( turn > 0 ) {
      GrapeDataDisplay< MyGridType > grape( gridPart.grid() );
      grape.dataDisplay( solution );
    }
  #endif

  interpolate( f, solution );
  double newL2error = l2norm.distance( f, solution );
  double newH1error = h1norm.distance( f, solution );

  std :: cout << "L2 error for interpolation after adaption: " << newL2error << std :: endl;
  std :: cout << "H1 error for interpolation after adaption: " << newH1error << std :: endl;

  #if USE_GRAPE && SHOW_INTERPOLATION
    if( turn > 0 ) {
      GrapeDataDisplay< MyGridType > grape( gridPart.grid );
      grape.dataDisplay( solution );
    }
  #endif

  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  std :: cout << "L2 EOC: " << l2eoc << std :: endl;
  std :: cout << "H1 EOC: " << h1eoc << std :: endl;

  /*
  const bool isLocallyAdaptive = Dune::Fem::Capabilities::isLocallyAdaptive< GridPartType :: GridType > :: v ;
  // threshold for EOC difference to predicted value
  const double eocThreshold = Parameter :: getValue("adapt.eocthreshold", double(0.2) );

  if( isLocallyAdaptive )
  {
    const double sign = step / std::abs( step );
    if( std::abs( l2eoc - h1eoc - sign ) > eocThreshold )
      DUNE_THROW( InvalidStateException,"Wrong L2/H1 relation");

    if( std::abs( l2eoc - ( sign * ( solution.space().order() + 1) ) ) > eocThreshold )
      DUNE_THROW( InvalidStateException,"Wrong L2-EOC for " << (step > 0) ? "refinement" : "coarsening" );
    if( std::abs( h1eoc - ( sign * ( solution.space().order() ) ) ) > eocThreshold )
      DUNE_THROW( InvalidStateException,"Wrong H1-EOC for " << (step > 0) ? "refinement" : "coarsening" );
  }
  */

  #if WRITE_DATA
    GrapeDataIO< MyGridType > dataio;
    dataio.writeGrid( gridPart.grid(), xdr, "gridout", 0, turn );
    dataio.writeData( solution, xdr, "sol", turn );
  #endif

  std :: cout << std :: endl;
}



int Main ( int argc, char **argv )
try
{
  MPIManager :: initialize( argc, argv );

  const char* paramName = "parameter";
  if( argc < 2 )
  {
    std :: cerr << "Usage: " << argv[ 0 ] << "<parameter>" << std :: endl;
  }
  else
    paramName = argv[1];

  std::string paramFile( paramName );

  // append parameter
  Parameter :: append( argc , argv );
  Parameter :: append( paramFile );

  int ml = 2 ; // default value = 2
  ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  MyGridType &grid = Dune::Fem::TestGrid::grid();
  const int step = Dune::Fem::TestGrid::refineStepsForHalf();

  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionSpaceTwoType discreteFunctionSpaceTwo( gridPart );
  DiscreteFunctionType solutionOne( "solution 1", discreteFunctionSpace );
  DiscreteFunctionTwoType solutionTwo( "solution 2", discreteFunctionSpaceTwo );
  DiscreteFunctionType solutionThree( "solution 3", discreteFunctionSpace );
  solutionOne.clear();
  solutionTwo.clear();
  solutionThree.clear();

  auto functionTuple = std::make_tuple(std::ref(solutionOne),
                                       std::ref(solutionTwo),
                                       std::ref(solutionThree));

  const bool locallyAdaptive = Parameter :: getValue< bool >("adapt.locallyadaptive", false );

  std :: cout << std :: endl << "Refining: " << std :: endl;
  for( int i = 0; i < ml; ++i )
    algorithm( gridPart, functionTuple, step, (i == ml-1), locallyAdaptive );

  std :: cout << std :: endl << "Coarsening:" << std::endl;
  for( int i = ml - 1; i >= 0; --i )
    algorithm( gridPart, functionTuple, -step, 1, locallyAdaptive );

  return 0;
}
catch( const Dune :: Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
