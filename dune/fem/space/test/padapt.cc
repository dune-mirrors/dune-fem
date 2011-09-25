#include <config.h>

#define SHOW_INTERPOLATION 1
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

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/operator/lagrangeinterpolation.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

// include solvers
#include <dune/fem/solver/inverseoperators.hh>
#include <dune/fem/solver/oemsolver.hh>

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

#include "systemmatrix.hh"


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
// typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//! type of the discrete function space our unkown belongs to
typedef PAdaptiveLagrangeSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;

//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef ExactSolution< FunctionSpaceType > ExactSolutionType;

typedef GridFunctionAdapter< ExactSolutionType, GridPartType >
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


void setPolOrder( const DiscreteFunctionSpaceType &space, bool increase ) 
{
  const int maxPol = 2;
  const int minPol = 1;
  const int size = space.indexSet().size( 0 );

  const int p = increase ? minPol : maxPol ;
  std::vector< int > polOrds( size, p );

  for( int i=0; i<size/2; ++i ) 
  {
    const int p = increase ? maxPol : minPol ;
    polOrds[ i ] = p;
  }

  std::cout << "Set polynomial order to " << p << std::endl;
  space.adapt( polOrds );
}

 

void adapt ( MyGridType &grid, DiscreteFunctionType &solution, int step )
{
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType;
  
  const DiscreteFunctionSpaceType &discreteFunctionSpace = solution.space();
  
  /*
  RestrictProlongOperatorType rp( solution );

  AdaptationManagerType adaptationManager( grid, rp );

  std :: string message = (step < 0 ? "Coarsening..." : "Refining..." );
  const int mark = (step < 0 ? -1 : 1);
  const int count = std :: abs( step );

  for( int i = 0; i < count; ++i ) 
  {
    IteratorType it = discreteFunctionSpace.begin();
    const IteratorType endit = discreteFunctionSpace.end();
    for( ; it != endit; ++it )
    {
      grid.mark( mark, *it );
    }

    // adapt grid 
    adaptationManager.adapt();

  }
  */

  setPolOrder( discreteFunctionSpace, step > 0 );
}

bool checkContinuous( DiscreteFunctionType &solution )
{
  double ret = 0;
  typedef DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef DiscreteFunctionSpaceType :: IteratorType IteratorType ;
  typedef IteratorType :: Entity  EntityType ;
  typedef GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef GridPartType :: IntersectionType         IntersectionType;
  
  const IteratorType endit = solution.space().end();
  for( IteratorType it = solution.space().begin(); it != endit; ++it ) 
  {
    const EntityType& entity = *it;
    const IntersectionIteratorType endiit = solution.space().gridPart().iend( entity );
    for( IntersectionIteratorType iit = solution.space().gridPart().ibegin( entity ); iit != endiit ; ++ iit ) 
    {
      const IntersectionType& intersection = *iit ; 
	    if( intersection.neighbor() && intersection.conforming() )
      {
        typedef CachingQuadrature< GridPartType, 1 > FaceQuadratureType;
        typedef IntersectionQuadrature< FaceQuadratureType, true > IntersectionQuadratureType;
        typedef IntersectionQuadratureType :: FaceQuadratureType QuadratureImp;
        IntersectionQuadratureType interQuad( solution.space().gridPart(), intersection, 4 );
        const QuadratureImp &quadInside  = interQuad.inside();
        const QuadratureImp &quadOutside = interQuad.outside();
        for( unsigned int qp = 0; qp < quadInside.nop(); ++qp )
	      {
          DiscreteFunctionType::RangeType uIn,uOut;
          solution.localFunction(entity).evaluate(quadInside[qp], uIn);
          solution.localFunction(*(intersection.outside())).evaluate(quadOutside[qp], uOut);
          ret = std::max(ret, (uIn-uOut).two_norm());
        }
      }
    }
  }
  return (ret<1e-10);
}

template <class Function>
void interpolate( const Function &f, DiscreteFunctionType &solution )
{
#if 0
  std::cout << "Applying Lagrangeinterpolation:" << std::endl;
  LagrangeInterpolation< DiscreteFunctionType > :: interpolateFunction( f, solution );
#else
  std::cout << "Applying L2-projection:" << std::endl;
  // define Laplace operator
  typedef MassModel< FunctionSpaceType > ModelType;
  typedef EllipticOperator< DiscreteFunctionType, ModelType > EllipticOperatorType;
  typedef Dune::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::GridPartType GridPartType;

  // typedef Dune::GridFunctionAdapter< Function, GridPartType > GridRHSFunctionType;
  // GridRHSFunctionType gridf( "grid rhs", f, solution.space().gridPart(), solution.space().order()+1 );
  DiscreteFunctionType rhs( "rhs", solution.space() );
  assembleRHS( f, rhs );
  
  EllipticOperatorType ellipticOp;

  // create linear inverse operator
  const double solverEps = Dune::Parameter::getValue< double >( "solvereps", 1e-8 );
  LinearInverseOperatorType solver( ellipticOp, solverEps, solverEps );

  solver( rhs, solution );
#endif
  bool continuous = checkContinuous( solution );

  L2Norm< GridPartType > l2norm( solution.space().gridPart() );
  H1Norm< GridPartType > h1norm( solution.space().gridPart() );
  double preL2error = l2norm.distance( f, solution );
  double preH1error = h1norm.distance( f, solution );

  std::cout << "Solution is " << (continuous?"":"NOT") << " continuous" << std::endl;
  std::cout << "L2 error before adaptation: " << preL2error << std::endl;
  std::cout << "H1 error before adaptation: " << preH1error << std::endl; 
}

void algorithm ( GridPartType &gridPart,
                 DiscreteFunctionType &solution, 
                 int step,
                 int turn )
{
  std::cout << "********************************************************" << std::endl;

  const unsigned int polOrder
    = DiscreteFunctionSpaceType :: polynomialOrder + 1;

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  interpolate( f, solution );
  
  std::cout << "Unknowns before adaptation: " << solution.space().size() << std::endl;
  adapt( gridPart.grid(), solution, step );
  
  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );
  double postL2error = l2norm.distance( f, solution );
  double postH1error = h1norm.distance( f, solution );

  //solution.print( std::cout );
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
    //if( turn > 0 ) 
    {
      GrapeDataDisplay< MyGridType > grape( gridPart.grid() );
      grape.dataDisplay( solution );
    }
  #endif

  interpolate( f, solution );
  
  #if USE_GRAPE && SHOW_INTERPOLATION
    //if( turn > 0 ) 
    {
      GrapeDataDisplay< MyGridType > grape( gridPart.grid() );
      grape.dataDisplay( solution );
    }
  #endif

/*
  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  std :: cout << "L2 EOC: " << l2eoc << std :: endl;
  std :: cout << "H1 EOC: " << h1eoc << std :: endl;
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
  //ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::ostringstream gridName;
  gridName << MyGridType::dimensionworld << "dgrid.dgf";
  GridPtr< MyGridType > gridptr( gridName.str().c_str() );

  const int step = 2;

  GridPartType gridPart( *gridptr );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  setPolOrder( discreteFunctionSpace, false );

  for ( int r = 0; r < 4; ++r )
  {
    std :: cout << std :: endl << "Refining: " << std :: endl;
    for( int i = 0; i < ml; ++i )
      algorithm( gridPart, solution, step, (i == ml-1) );
    
    std :: cout << std :: endl << "Coarsening:" << std::endl;
    for( int i = ml - 1; i >= 0; --i )
      algorithm( gridPart, solution, -step, 1 );

    if (r==0)
    {
      // Test grid ref. (polynomial order is max)
      std :: cout << std :: endl << "Refine grid" << std::endl;
      Dune::GlobalRefine::apply(*gridptr,1);
    }
		else if (r==1)
    {
      // Test grid ref. (but with polynomial order set to min)
      std :: cout << std :: endl << "Refine grid" << std::endl;
      setPolOrder( discreteFunctionSpace, true );
      Dune::GlobalRefine::apply(*gridptr,1);
    }
    else if (r==2)
    {
      // Test grid coarsening (polynomial order set to max)
      std :: cout << std :: endl << "Coarsen grid" << std::endl;
      // Dune::GlobalRefine::apply(*gridptr,-1); // not workling yet
    }
  }

  return 0;
}
catch( const Dune :: Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
