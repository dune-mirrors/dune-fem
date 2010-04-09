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
#include <dune/common/stdstreams.cc>

#include <dune/fem/gridpart/gridpart.hh>
#include <dgfgridtype.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/space/lagrangespace/adaptmanager.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
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
#include <dune/fem/io/parameter.hh>


// Check for unhealthy grids
// -------------------------

// forward declaration of the real main method
int Main ( int argc, char **argv );

template< class Grid >
struct CheckGridEnabled
{
  typedef Grid GridType;

  //typedef Dune::HierarchicGridPart< GridType > GridPartType;
  typedef Dune::AdaptiveLeafGridPart< GridType > GridPartType;
  
  inline static int CallMain ( int argc, char **argv )
  {
    return Main( argc, argv );
  }
};

// disable YaspGrid
namespace Dune
{
#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
  template< int dim >
  class YaspGrid;
#else
  template< int dim, int dimworld >
  class YaspGrid;
#endif
}

#if DUNE_VERSION_NEWER(DUNE_GRID,1,3,0)
template< int dim >
struct CheckGridEnabled< Dune :: YaspGrid< dim > >
{
  typedef Dune :: YaspGrid< dim > GridType;
  
  typedef Dune :: LeafGridPart< GridType > GridPartType;

  inline static int CallMain ( int argc, char **argv )
  {
    std :: cerr << "WARNING: Lagrange Adaptation test disabled, because YaspGrid sucks!"
                << std :: endl;
    return 0;
  }
};
#else
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
#endif

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
  return CheckGridEnabled< Dune::GridSelector::GridType >::CallMain( argc, argv );
}



using namespace Dune;

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

//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef DiscreteFunctionAdapter< ExactSolutionType, GridPartType >
  GridExactSolutionType;

//! type of the DoF manager
typedef DofManager< MyGridType > DofManagerType;
//! type of the DoF manager factory
typedef DofManagerFactory< DofManagerType > DofManagerFactoryType;

//! type of restrict-prolong operator
typedef RestrictProlongDefault< DiscreteFunctionType >
  RestrictProlongOperatorType;
//! type of the adaption manager
typedef AdaptationManager< MyGridType, RestrictProlongOperatorType >
  AdaptationManagerType;




 

void adapt ( MyGridType &grid, DiscreteFunctionType &solution, int step )
{
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  
  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();
  
  RestrictProlongOperatorType rp( solution );

  AdaptationManagerType adaptationManager( grid, rp );

  std :: string message = (step < 0 ? "Coarsening..." : "Refining..." );
  const int mark = (step < 0 ? -1 : 1);
  const int count = std :: abs( step );
  
  for( int i = 0; i < count; ++i ) {
    IteratorType it = discreteFunctionSpace.begin();
    const IteratorType endit = discreteFunctionSpace.end();
    for( ; it != endit; ++it )
      grid.mark( mark, *it );

    // adapt grid 
    adaptationManager.adapt();
  }
}



void algorithm ( GridPartType &gridPart,
                 DiscreteFunctionType &solution, 
                 int step,
                 int turn )
{
  const unsigned int polOrder
    = DiscreteFunctionSpaceType :: polynomialOrder + 1;

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );
  
  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
  double preL2error = l2norm.distance( f, solution );
  double preH1error = h1norm.distance( f, solution );

  std::cout << "Unknowns before adaptation: " << solution.space().size() << std::endl;
  std::cout << "L2 error before adaptation: " << preL2error << std::endl;
  std::cout << "H1 error before adaptation: " << preH1error << std::endl; 
  
  adapt( gridPart.grid(), solution, step );
  
  double postL2error = l2norm.distance( f, solution );
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

  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
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

  std::ostringstream gridName;
  gridName << MyGridType::dimensionworld << "dgrid.dgf";
  GridPtr< MyGridType > gridptr( gridName.str().c_str() );

  const int step = DGFGridInfo< MyGridType >::refineStepsForHalf();

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
catch( const Dune :: Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
