#include <config.h>

// to write out the data, set WRITE_DATA to 1
#define WRITE_DATA 0

#ifndef USE_LFE
#define USE_LFE 0
#endif

#define USE_PSPACE 0

// polynomial order of base functions
const int polOrder = 2; // POLORDER;

#include <iostream>
#include <sstream>
//#include <dune/common/stdstreams.cc>

#include <dune/fem/gridpart/leafgridpart.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/io/file/dataoutput.hh>

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

template< int dim, class CoordCont >
struct CheckGridEnabled< Dune :: YaspGrid< dim, CoordCont > >
{
  typedef Dune :: YaspGrid< dim, CoordCont > GridType;

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
    {
      phi[ 0 ] *= sin( M_PI * x[ i ] );
      phi[ 1 ] += x[ i ] * x[ i ];
    }
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
      {
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
        Dphi[ 1 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
      }
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
typedef FunctionSpace< double, double, MyGridType::dimensionworld, 2 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
#if USE_LFE
#warning "Using LocalFiniteElement based LagrangeSpace!"
typedef LagrangeSpace< FunctionSpaceType, GridPartType >
  DiscreteFunctionSpaceType;
#elif USE_PSPACE
typedef Fem :: PAdaptiveLagrangeSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;
#else
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;
#endif

//! type of the discrete function we are using
#if HAVE_PETSC
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#else
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
#endif

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef GridFunctionAdapter< ExactSolutionType, GridPartType >
  GridExactSolutionType;

//! type of the DoF manager
typedef DofManager< MyGridType > DofManagerType;

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< MyGridType, RestrictProlongOperatorType >
  AdaptationManagerType;






void adapt ( MyGridType &grid, DiscreteFunctionType &solution, int step,
             const bool locallyAdaptive )
{
  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();

  RestrictProlongOperatorType rp( solution );

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



void algorithm ( GridPartType &gridPart,
                 DiscreteFunctionType &solution,
                 int step,
                 int turn,
                 const bool locallyAdaptive )
{
  const unsigned int polOrder = solution.space().order() + 1;

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder+1 );

  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );

  interpolate( f, solution );
  double preL2error = l2norm.distance( f, solution );
  double preH1error = h1norm.distance( f, solution );

  size_t unknowns = gridPart.comm().sum( solution.space().size() );
  if( Parameter::verbose() )
  {
    std::cout << "Unknowns before adaptation: " << unknowns   << std::endl;
    std::cout << "L2 error before adaptation: " << preL2error << std::endl;
    std::cout << "H1 error before adaptation: " << preH1error << std::endl;
  }

  adapt( gridPart.grid(), solution, step, locallyAdaptive );

  // check index set for being consistent with the dune index set check
  checkAdaptiveIndexSet( gridPart );

  //double postL2error = l2norm.distance( fexact, solution );
  double postL2error = l2norm.distance( solution, fexact );
  double postH1error = h1norm.distance( f, solution );

  unknowns = gridPart.comm().sum( solution.space().size() );
  if( Parameter::verbose() )
  {
    std::cout << "Unknowns after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << unknowns << std::endl;
    std::cout << "L2 error after "
                << (step < 0 ? "restriction" : "prolongation")
                << ": " << postL2error << std::endl;
    std::cout << "H1 error after "
                << (step < 0 ? "restriction" : "prolongation")
                << ": " << postH1error << std::endl;

  }
#if WRITE_DATA
  {
    static int turn = 0;
    typedef std::tuple< DiscreteFunctionType* > IODataType;
    IODataType data( &solution );
    Dune::Fem::DataOutput< MyGridType, IODataType > output( gridPart.grid(), data );
    output.writeData( turn, "test" );
    ++turn;
  }
#endif

  interpolate( f, solution );
  double newL2error = l2norm.distance( f, solution );
  double newH1error = h1norm.distance( f, solution );

  if( Parameter::verbose() )
  {
    std :: cout << "L2 error for interpolation after adaption: " << newL2error << std :: endl;
    std :: cout << "H1 error for interpolation after adaption: " << newH1error << std :: endl;
  }

  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  if( Parameter::verbose() )
  {
    std :: cout << "L2 EOC: " << l2eoc << std :: endl;
    std :: cout << "H1 EOC: " << h1eoc << std :: endl;
  }

  const bool isLocallyAdaptive = Dune::Fem::Capabilities::isLocallyAdaptive< GridPartType :: GridType > :: v ;
  // threshold for EOC difference to predicted value
  const double eocThreshold = Parameter :: getValue("adapt.eocthreshold", double(0.25) );

  if( isLocallyAdaptive )
  {
    const double sign = step / std::abs( step );
    if( std::abs( l2eoc - h1eoc - sign ) > eocThreshold )
      DUNE_THROW( InvalidStateException,"Wrong L2/H1 relation");

    const char* refcrs = (step > 0) ? "refinement" : "coarsening";
    if( std::abs( l2eoc - ( sign * ( solution.space().order() + eocThreshold) ) ) < 0 )
      DUNE_THROW( InvalidStateException,"Wrong L2-EOC for " << refcrs );
    if( std::abs( h1eoc - ( sign * ( solution.space().order()-1.0+eocThreshold ) ) ) < 0 )
      DUNE_THROW( InvalidStateException,"Wrong H1-EOC for " << refcrs );
  }


  if( Parameter::verbose() )
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
#if USE_LFE
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart, polOrder );
#else
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
#endif
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  const bool locallyAdaptive = Parameter :: getValue< bool >("adapt.locallyadaptive", false );

  if( Parameter::verbose() )
    std :: cout << std :: endl << "Refining: " << std :: endl;
  for( int i = 0; i < ml; ++i )
    algorithm( gridPart, solution, step, (i == ml-1), locallyAdaptive );

  if( Parameter::verbose() )
    std :: cout << std :: endl << "Coarsening:" << std::endl;
  for( int i = ml - 1; i >= 0; --i )
    algorithm( gridPart, solution, -step, 1, locallyAdaptive );

  return 0;
}
catch( const Dune :: Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
