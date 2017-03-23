#include <config.h>

// C++ includes
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-grid includes
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/fourier.hh>
#include <dune/fem/space/lagrange.hh>

#include "../../test/exactsolution.hh"


// polynomial order
static const int polOrder = POLORDER;

// range dimension
static const int dimRange = DIMRANGE;

// grid type
typedef Dune::GridSelector::GridType GridType;



// interpolationError
// ------------------

template< class DiscreteFunctionSpace >
double interpolationError ( const DiscreteFunctionSpace &discreteFunctionSpace )
{
  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
  const GridPartType &gridPart = discreteFunctionSpace.gridPart();

  // create exact solution
  typedef Dune::Fem::ExactSolution< typename DiscreteFunctionSpace::FunctionSpaceType > ExactSolutionType;
  ExactSolutionType exactSolution;
  typedef Dune::Fem::GridFunctionAdapter< ExactSolutionType, GridPartType > GridExactSolutionType;
  GridExactSolutionType gridExactSolution( "exact solution", exactSolution, gridPart, 5 );

  // create discrete function
  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  // perform the interpolation
  interpolate( gridExactSolution, solution );

  typedef Dune::Capabilities::hasSingleGeometryType< typename GridPartType::GridType > HasSingleGeometryType;
  if( (DiscreteFunctionSpace::dimDomain == 2) && HasSingleGeometryType::v && (HasSingleGeometryType::topologyId == 0) )
  {
    typename DiscreteFunctionType::HessianRangeType h1, h2;
    typename DiscreteFunctionType::DomainType x( 0.5 );

    solution.hessian( x, h1 );
    exactSolution.hessian( x, h2 );

    double norm =0;
    for( int r =0; r< dimRange; ++r )
    {
      h1[ r ]-= h2[ r ];
      norm += h1[r ].frobenius_norm2();
    }

    std::cout<<std::sqrt( norm )<<std::endl;
  }

  // compute error
  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  return l2norm.distance( gridExactSolution, solution );
}



// interpolationErrors
// -------------------

template< class... DiscreteFunctionSpace >
std::array< double, sizeof...( DiscreteFunctionSpace ) > interpolationErrors ( const DiscreteFunctionSpace &... discreteFunctionSpace )
{
  return {{ interpolationError( discreteFunctionSpace )... }};
}



// printValues
// -----------

template< std::size_t n >
void printValues ( std::string prefix, std::array< double, n > values )
{
  if( Dune::Fem::MPIManager::rank() != 0 )
    return;

  std::cout << prefix << "( " << values[ 0 ];
  for( std::size_t i = 1; i < n; ++i )
    std::cout << ", " << values[ i ];
  std::cout << " )" << std::endl;
}



// eocLoop
// -------

template< class... DiscreteFunctionSpace >
void eocLoop ( GridType &grid, int steps, const DiscreteFunctionSpace &... discreteFunctionSpace )
{
  const std::size_t n = sizeof...( DiscreteFunctionSpace );
  std::array< double, n > oldErrors = interpolationErrors( discreteFunctionSpace... );
  printValues( "L2 error[ 0 ] = ", oldErrors );

  for( int step = 1; step <= steps ; ++step )
  {
    Dune::Fem::GlobalRefine::apply( grid, Dune::DGFGridInfo< GridType >::refineStepsForHalf() );

    const std::array< double, n > newErrors = interpolationErrors( discreteFunctionSpace... );
    std::array< double, n > eocs;
    std::transform( oldErrors.begin(), oldErrors.end(), newErrors.begin(), eocs.begin(), [] ( double eold, double enew ) { return log( eold / enew ) / M_LN2; } );
    oldErrors = newErrors;

    printValues( "L2 error[ " + std::to_string( step ) + " ] = ", oldErrors );
    printValues( "L2 EOC[ " + std::to_string( step ) + " ] = ", eocs );
  }
}



// main
// ----

int main ( int argc, char **argv )
try
{
  // initialize MPI
  Dune::Fem::MPIManager::initialize( argc, argv );

  // read parameter file
  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  if( argc < 1 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " [initial refinement level] [steps]" << std::endl;
    return 0;
  }

  // initial refinement level
  const int refCount = (argc > 1 ? std::max( 0, std::atoi( argv[ 1 ] ) ) : 1);

  // number of EOC steps
  const int steps = (argc > 2 ? std::max( 0, std::atoi( argv[ 2 ] ) ) : 1);

  // create grid
  Dune::GridPtr< GridType > grid( std::to_string( GridType::dimension ) + "dgrid.dgf" );
  grid->loadBalance();
  Dune::Fem::GlobalRefine::apply( *grid, refCount * Dune::DGFGridInfo< GridType >::refineStepsForHalf() );

  // create grid part
  typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  // function space type
  typedef typename GridPartType::ctype DomainFieldType;
  static const int dimDomain = GridPartType::dimensionworld;
  typedef Dune::Fem::FunctionSpace< DomainFieldType, double, dimDomain, dimRange > FunctionSpaceType;

  // create discrete function spaces
  Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder > discontinuousGalerkinSpace( gridPart );
  Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > lagrangeSpace( gridPart );
  Dune::Fem::FourierDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 > fourierSpace( gridPart, polOrder+1 );

  // perform eoc loop
  eocLoop( *grid, steps, discontinuousGalerkinSpace, lagrangeSpace, fourierSpace );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
