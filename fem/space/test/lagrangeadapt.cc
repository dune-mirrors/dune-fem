#include <config.h>

#define SHOW_INTERPOLATION 0
#define SHOW_RESTRICT_PROLONG 1

// to write out the data, set WRITE_DATA to 1
#define WRITE_DATA 0

// to use generic adaption, set GENERIC_ADAPT to 1
#define GENERIC_ADAPT 1

// to use grape, set to WANT_GRAPE to 1
#ifndef WANT_GRAPE
#define WANT_GRAPE 0
#endif

// polynomial order of base functions
const int polOrder = POLORDER;

#include <iostream>
#include <dune/common/stdstreams.cc>

#include <dune/fem/gridpart/gridpart.hh>
#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include <dune/fem/space/common/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/adaptmanager.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachequad.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
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
#include <dune/fem/io/file/grapedataio.hh>


// Check for unhealthy grids
// -------------------------

// forward declaration of the real main method
int Main ( int argc, char **argv );

template< class Grid >
struct CheckGridEnabled
{
  typedef Grid GridType;

  typedef Dune :: AdaptiveLeafGridPart< GridType > GridPartType;
  
  inline static int CallMain ( int argc, char **argv )
  {
    return Main( argc, argv );
  }
};

// disable YaspGrid
namespace Dune
{
  template< int dim, int dimworld >
  class YaspGrid;
}

template< int dim, int dimworld >
struct CheckGridEnabled< Dune :: YaspGrid< dim, dimworld > >
{
  typedef Dune :: YaspGrid< dim, dimworld > GridType;
  
  typedef Dune :: LeafGridPart< GridType > GridPartType;

  inline static int CallMain ( int argc, char **argv )
  {
    std :: cerr << "WARNING: Lagrange Adaptation test disabled, because YaspGrid sucks!"
                << std :: endl;
    return 0;
  }
};

// disable UGGrid
namespace Dune
{
  template< int dim >
  class UGGrid;
}

template< int dim >
struct CheckGridEnabled< Dune :: UGGrid< dim > >
{
  typedef Dune :: UGGrid< dim > GridType;

  typedef Dune :: LeafGridPart< GridType > GridPartType;

  inline static int CallMain ( int argc, char **argv )
  {
    std :: cerr << "WARNING: Lagrange Adaptation test disabled, because UGGrid sucks!"
                << std :: endl;
    return 0;
  }
};

int main ( int argc, char **argv )
{
  return CheckGridEnabled< GridType > :: CallMain( argc, argv );
}



using namespace Dune;

// Function, we will interpolate
// -----------------------------

template< class FunctionSpace >
class ExactSolution
: public Function< FunctionSpace, ExactSolution< FunctionSpace > >
{
public:
  typedef FunctionSpace FunctionSpaceType;

private:
  typedef ExactSolution< FunctionSpaceType > ThisType;
  typedef Function< FunctionSpaceType, ThisType > BaseType;

public:
  typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
  typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

  typedef typename FunctionSpaceType :: DomainType DomainType;
  typedef typename FunctionSpaceType :: RangeType RangeType;
  typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;

public:
  ExactSolution ( FunctionSpaceType &functionSpace )
  : BaseType( functionSpace )
  {
  }

  void evaluate ( const DomainType &x,
                  RangeType &phi ) const
  {
    phi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      // phi[ 0 ] += x[ i ] * x[ i ]; 
      phi[ 0 ] *= sin( M_PI * x[ i ] ); 
  }

  void evaluate ( const DomainType &x,
                  RangeFieldType t,
                  RangeType &phi ) const
  {
    evaluate( x, phi );
  }

  void jacobian( const DomainType &x,
                 JacobianRangeType &Dphi ) const
  {
    Dphi = 1;
    for( int i = 0; i < DomainType :: dimension; ++i )
      for( int j = 0; j < DomainType :: dimension; ++j )
        // Dphi[ 0 ][ j ] *= ((i != j) ? 1. : 2.*x[i]);
        Dphi[ 0 ][ j ] *= ((i != j) ? sin( M_PI * x[ i ]) : M_PI * cos( M_PI * x[ i ] ));
  }

  void jacobian( const DomainType &x,
                 RangeFieldType t,
                 JacobianRangeType &Dphi ) const
  {
    jacobian( x, Dphi );
  }
};



// Type Definitions
// ----------------

typedef CheckGridEnabled< GridType > :: GridPartType GridPartType;

//! type of the function space
typedef FunctionSpace< double, double, dimworld, 1 > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
  DiscreteFunctionSpaceType;

//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
  GridExactSolutionType;

//! type of the DoF manager
typedef DofManager< GridType > DofManagerType;
//! type of the DoF manager factory
typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< GridType, RestrictProlongOperatorType >
  AdaptationManagerType;




 

void adapt ( GridType &grid, DiscreteFunctionType &solution, int step )
{
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  
  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();
  
  RestrictProlongOperatorType rp( solution );

  #if GENERIC_ADAPT 
    AdaptationManagerType adaptationManager( grid, rp );
  #else 
    DofManagerType &dofManager = DofManagerFactoryType :: getDofManager( grid );
  #endif  

  std :: string message = (step < 0 ? "Coarsening..." : "Refining..." );
  const int mark = (step < 0 ? -1 : 1);
  const int count = std :: abs( step );
  
  for( int i = 0; i < count; ++i ) {
    IteratorType it = discreteFunctionSpace.begin();
    const IteratorType endit = discreteFunctionSpace.end();
    for( ; it != endit; ++it )
      grid.mark( mark, it );
      
    #if GENERIC_ADAPT
      std :: cout << message << std::endl;
      adaptationManager.adapt();
    #else
      std :: cout << message << " (nongeneric)" << std :: endl;
      grid.adapt( dofManager, rp );
    #endif
  }
}



void algorithm ( GridPartType &gridPart,
                 DiscreteFunctionType &solution, 
                 int step,
                 int turn )
{
  const unsigned int polOrder
    = DiscreteFunctionSpaceType :: polynomialOrder + 1;
  
  FunctionSpaceType functionSpace;
  ExactSolutionType fexact( functionSpace );
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );
  
  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
  double preL2error = l2norm.distance( f, solution );
  double preH1error = h1norm.distance( f, solution );

  std :: cout << "L2 error before adaption: " << preL2error << std :: endl;
  std :: cout << "H1 error before adaption: " << preH1error << std :: endl; 
  
  adapt( gridPart.grid(), solution, step );
  
  double postL2error = l2norm.distance( f, solution );
  double postH1error = h1norm.distance( f, solution );

  std :: cout << "L2 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postL2error << std :: endl;
  std :: cout << "H1 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postH1error << std :: endl; 
  
  #if USE_GRAPE && SHOW_RESTRICT_PROLONG
    if( turn > 0 ) {
      GrapeDataDisplay< GridType > grape( gridPart.grid() );
      grape.dataDisplay( solution );
    }
  #endif

  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
  double newL2error = l2norm.distance( f, solution );
  double newH1error = h1norm.distance( f, solution );

  std :: cout << "L2 error for interpolation after adaption: " << newL2error << std :: endl;
  std :: cout << "H1 error for interpolation after adaption: " << newH1error << std :: endl; 
  
  #if USE_GRAPE && SHOW_INTERPOLATION
    if( turn > 0 ) {
      GrapeDataDisplay< GridType > grape( gridPart.grid );
      grape.dataDisplay( solution );
    }
  #endif

  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  std :: cout << "L2 EOC: " << l2eoc << std :: endl;
  std :: cout << "H1 EOC: " << h1eoc << std :: endl;

  #if WRITE_DATA
    GrapeDataIO< GridType > dataio; 
    dataio.writeGrid( gridPart.grid(), xdr, "gridout", 0, turn );
    dataio.writeData( solution, xdr, "sol", turn );
  #endif

  std :: cout << std :: endl;
}



int Main ( int argc, char **argv )
{
  if( argc != 2 ) {
    std :: cerr << "Usage: " << argv[ 0 ] << "<maxlevel>" << std :: endl;
    exit( 1 );
  }
  
  int ml = atoi( argv[ 1 ] );

  char tmp[ 100 ]; 
  sprintf( tmp, "%ddgrid.dgf", dimworld );
  GridPtr< GridType > gridptr( tmp );

  const int step = DGFGridInfo< GridType > :: refineStepsForHalf();

  GridPartType gridPart( *gridptr );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  std :: cout << std :: endl << "Refining: " << std :: endl;
  for( int i = 0; i < ml; ++i )
    algorithm( gridPart, solution, step, (i == ml-1) );
  
  std :: cout << std :: endl << "Coarsening:" << std::endl;
  for( int i = ml - 1; i >= 0; --i )
    algorithm( gridPart, solution, -step, 1 );

  return 0;
}
