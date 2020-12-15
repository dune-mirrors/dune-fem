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
//#include <dune/common/stdstreams.cc>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>

#include <dune/fem/quadrature/intersectionquadrature.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/cachingquadrature.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/common/dataprojection/tuple.hh>
#include <dune/fem/space/hpdg/orthogonal.hh>

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

int main ( int argc, char **argv )
{
  return CheckGridEnabled< Dune::GridSelector::GridType >::CallMain( argc, argv );
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
typedef FunctionSpace< double, double, MyGridType::dimensionworld, MyGridType::dimensionworld > FunctionSpaceType;

//! type of the discrete function space our unkown belongs to
// typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//! type of the discrete function space our unkown belongs to
//typedef LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
typedef Fem :: PAdaptiveLagrangeSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
//typedef Fem :: PAdaptiveDGSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;

//typedef hpDG::OrthogonalDiscontinuousGalerkinSpace< FunctionSpaceType,
//                                    GridPartType, polOrder, true >  DiscreteFunctionSpaceType;



//! type of the discrete function we are using
typedef AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
typedef Dune::Fem::DefaultDataProjectionTuple< DiscreteFunctionType > DataProjectionType;

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


void setPolOrder( DiscreteFunctionSpaceType &space, DiscreteFunctionType& solution, bool increase )
{
  const int maxPol = POLORDER;
  const int minPol = 1;
  const int size = space.indexSet().size( 0 );

  const int p = increase ? minPol : maxPol ;
  std::vector< int > polOrds( size, p );

  for( int i=0; i<size/2; ++i )
  {
    const int p = increase ? maxPol : minPol ;
    polOrds[ i ] = p;
  }

  // mark space for all entities
  const auto end = space.end();
  for( auto it = space.begin(); it != end ; ++it )
  {
    const auto& entity = *it;
    space.mark( polOrds[ space.indexSet().index( entity )], entity );
  }

  Dune::Fem::hpDG::AdaptationManager< DiscreteFunctionSpaceType,
    DataProjectionType > pAdaptManager( const_cast< DiscreteFunctionSpaceType& >(space), DataProjectionType( solution ) );

  std::cout << "Set polynomial order to " << p << std::endl;
  pAdaptManager.adapt();
}

void polOrderAdapt( DiscreteFunctionSpaceType& space, DiscreteFunctionType &solution, int step)
{
  setPolOrder( space, solution, step > 0 );
}

void gridAdapt( MyGridType &grid, DiscreteFunctionType &solution, int step,
                const bool locallyAdaptive = false )
{
  #if USE_GRAPE && SHOW_RESTRICT_PROLONG
    //if( turn > 0 )
    {
      GrapeDataDisplay< MyGridType > grape( grid );
      grape.dataDisplay( solution );
    }
  #endif
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

    int elemNo = 0;
    for( const auto& entity : discreteFunctionSpace)
    {
      if( elemNo < numElements )
        grid.mark( mark, entity );
      ++elemNo;
    }

    // adapt grid
    adaptationManager.adapt();
  }

  #if USE_GRAPE && SHOW_RESTRICT_PROLONG
    //if( turn > 0 )
    {
      GrapeDataDisplay< MyGridType > grape( grid );
      grape.dataDisplay( solution );
    }
  #endif

}

bool checkContinuous( DiscreteFunctionType &solution )
{
  double ret = 0;
  typedef DiscreteFunctionType :: DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef DiscreteFunctionSpaceType :: GridPartType GridPartType;
  typedef GridPartType :: IntersectionIteratorType IntersectionIteratorType;
  typedef GridPartType :: IntersectionType         IntersectionType;

  Dune::Fem::ConstLocalFunction< DiscreteFunctionType > uInside( solution );
  Dune::Fem::ConstLocalFunction< DiscreteFunctionType > uOutside( solution );

  for( const auto& entity : solution.space() )
  {
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

        const auto nb = intersection.outside();

        auto iGuard = bindGuard( uInside, entity);
        auto oGuard = bindGuard( uOutside, nb );

        for( unsigned int qp = 0; qp < quadInside.nop(); ++qp )
	      {
          typename DiscreteFunctionType::RangeType uIn,uOut;

          uInside.evaluate(quadInside[qp], uIn);
          uOutside.evaluate(quadOutside[qp], uOut);

          ret = std::max(ret, (uIn-uOut).two_norm());
        }
      }
    }
  }
  return (ret<1e-10);
}

template <class Function>
void interpolateSolution( const Function &f, DiscreteFunctionType &solution )
{
#if 1
  std::cout << "Applying Lagrangeinterpolation:" << std::endl;
  interpolate( f, solution );
#else
  std::cout << "Applying L2-projection:" << std::endl;
  // define Laplace operator
  typedef MassModel< FunctionSpaceType > ModelType;
  typedef EllipticOperator< DiscreteFunctionType, ModelType > EllipticOperatorType;
  typedef Dune::CGInverseOperator< DiscreteFunctionType > LinearInverseOperatorType;
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

  typedef Dune::SparseRowMatrixTraits < DiscreteFunctionSpaceType, DiscreteFunctionSpaceType > MatrixObjectTraits;
  typedef LagrangeMatrixTraits< MatrixObjectTraits > MatrixTraits;
  typedef Dune::SparseRowLinearOperator< DiscreteFunctionType, DiscreteFunctionType, MatrixTraits >
      LinearOperatorType;

  DiscreteFunctionType rhs( "rhs", solution.space() );
  assembleRHS( f, rhs );

  EllipticOperatorType ellipticOp;

  // create linear inverse operator
  const double solverEps = Dune::Parameter::getValue< double >( "solvereps", 1e-8 );

#if 1
  LinearInverseOperatorType solver( ellipticOp, solverEps, solverEps );
#else
  LinearOperatorType linearOp( "assempled elliptic operator", solution.space(), solution.space() );
  ellipticOp.jacobian( solution, linearOp );
  LinearInverseOperatorType solver( linearOp, solverEps, solverEps );
#endif

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
                 DiscreteFunctionSpaceType& space,
                 DiscreteFunctionType &solution,
                 int step,
                 int turn )
{
  std::cout << "********************************************************" << std::endl;

  const unsigned int polOrder
    = DiscreteFunctionSpaceType :: polynomialOrder + 1;

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  interpolateSolution( f, solution );

  std::cout << "Unknowns before adaptation: " << space.size() << std::endl;
  polOrderAdapt( space, solution, step );

  L2Norm< GridPartType > l2norm( gridPart );
  H1Norm< GridPartType > h1norm( gridPart );
  double postL2error = l2norm.distance( f, solution );
  double postH1error = h1norm.distance( f, solution );

  //solution.print( std::cout );
  std::cout << "Unknowns after "
            << (step < 0 ? "restriction" : "prolongation")
            << ": " << space.size() << std::endl;
  std::cout << "L2 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postL2error << std::endl;
  std::cout << "H1 error after "
              << (step < 0 ? "restriction" : "prolongation")
              << ": " << postH1error << std::endl;

  interpolateSolution( f, solution );

/*
  double l2eoc = -log( newL2error / preL2error) / M_LN2;
  double h1eoc = -log( newH1error / preH1error) / M_LN2;

  std :: cout << "L2 EOC: " << l2eoc << std :: endl;
  std :: cout << "H1 EOC: " << h1eoc << std :: endl;
*/


#if WRITE_DATA
  typedef std::tuple< DiscreteFunctionType* > IODataType;
  IODataType data( &solution );
  Dune::Fem::DataOutput< GridType, IODataType > output( grid, data );

  output.writeStep( turn );
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

  int ml = 2 ; // default value = 1
  ml = Parameter :: getValue ("lagrangeadapt.maxlevel", ml);

  std::stringstream gridFile;
  gridFile << "DGF" << std::endl;
  gridFile << "Interval" << std::endl;
  for( int d=0; d<MyGridType::dimensionworld; ++d )
    gridFile << "0 ";
  gridFile << std::endl;
  for( int d=0; d<MyGridType::dimensionworld; ++d )
    gridFile << "1 ";
  gridFile << std::endl;
  for( int d=0; d<MyGridType::dimensionworld; ++d )
    gridFile << "1 ";
  gridFile << std::endl;
  gridFile << "#" <<  std::endl;

  GridPtr< MyGridType > gridptr( gridFile );

  const int step = 2;

  GridPartType gridPart( *gridptr );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  setPolOrder( discreteFunctionSpace, solution, false );

  ExactSolutionType fexact;
  GridExactSolutionType f( "exact solution", fexact, gridPart, polOrder );

  const bool locallyAdaptive = Parameter :: getValue< bool >("adapt.locallyadaptive", false );

  interpolateSolution( f, solution );
  for ( int r = 0; r < step; ++r )
  {
    std :: cout << std :: endl << "Refining: " << std :: endl;
    for( int i = 0; i < ml; ++i )
      algorithm( gridPart, discreteFunctionSpace, solution, step, (i == ml-1) );

    std :: cout << std :: endl << "Coarsening:" << std::endl;
    for( int i = ml - 1; i >= 0; --i )
      algorithm( gridPart, discreteFunctionSpace, solution, -step, 1 );

    if (r==0)
    {
      // Test grid ref. (polynomial order is max)
      std :: cout << std :: endl << "Refine grid" << std::endl;
      gridAdapt( *gridptr, solution, step, locallyAdaptive ) ;
    }
		else if (r==1)
    {
      // Test grid ref. (but with polynomial order set to min)
      std :: cout << std :: endl << "Refine grid" << std::endl;
      setPolOrder( discreteFunctionSpace, solution, true );
      gridAdapt( *gridptr, solution, step, locallyAdaptive ) ;
    }
    else if (r==2)
    {
      // Test grid coarsening (polynomial order set to max)
      std :: cout << std :: endl << "Coarsen grid" << std::endl;
      gridAdapt( *gridptr, solution, -step, locallyAdaptive ) ;
    }
  }

  return 0;
}
catch( const Dune :: Exception &e )
{
  std :: cerr << e << std :: endl;
  return 1;
}
