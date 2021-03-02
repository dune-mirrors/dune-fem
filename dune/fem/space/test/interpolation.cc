#include <config.h>

// C++ includes
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>

// dune-grid includes
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

// dune-fem includes
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/fourier.hh>
#include <dune/fem/space/lagrange.hh>

#if HAVE_DUNE_LOCALFUNCTIONS
#include <dune/fem/space/localfiniteelement/quadratureinterpolation.hh>
#endif

#include "../../test/exactsolution.hh"


// polynomial order
static const int polOrder = POLORDER;

// range dimension
static const int dimRange = DIMRANGE;

// grid type
typedef Dune::GridSelector::GridType GridType;



// interpolationError
// ------------------

template< class ExactSolutionType, class DiscreteFunctionSpace >
double interpolationErrorImpl ( const ExactSolutionType& exactSolution, const DiscreteFunctionSpace &discreteFunctionSpace )
{
  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
  const GridPartType &gridPart = discreteFunctionSpace.gridPart();

  // create exact solution
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

template< class DiscreteFunctionSpace >
double interpolationError ( const DiscreteFunctionSpace &discreteFunctionSpace )
{
  // create exact solution
  typedef Dune::Fem::ExactSolution< typename DiscreteFunctionSpace::FunctionSpaceType > ExactSolutionType;
  ExactSolutionType exactSolution;
  return interpolationErrorImpl( exactSolution, discreteFunctionSpace );
}

template< class DiscreteFunctionSpace >
double interpolationErrorPolynomial ( const DiscreteFunctionSpace &discreteFunctionSpace )
{
  // create exact solution
  typedef Dune::Fem::Polynomial< typename DiscreteFunctionSpace::FunctionSpaceType > ExactSolutionType;
  ExactSolutionType polynomial( polOrder );
  return interpolationErrorImpl( polynomial, discreteFunctionSpace );
}


// interpolationErrors
// -------------------

template< class... DiscreteFunctionSpace >
std::array< double, sizeof...( DiscreteFunctionSpace ) > interpolationErrors ( const DiscreteFunctionSpace &... discreteFunctionSpace )
{
  return {{ interpolationError( discreteFunctionSpace )... }};
}

template< class... DiscreteFunctionSpace >
std::array< double, sizeof...( DiscreteFunctionSpace ) > interpolationErrorsPolynomial ( const DiscreteFunctionSpace &... discreteFunctionSpace )
{
  return {{ interpolationErrorPolynomial( discreteFunctionSpace )... }};
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
  std::array< double, n > oldErrors  = interpolationErrors( discreteFunctionSpace... );
  std::array< double, n > oldErrors2 = interpolationErrorsPolynomial( discreteFunctionSpace... );
  printValues( "L2 error[ 0 ] = ", oldErrors );
  printValues( "L2 error[ 0 ] = ", oldErrors2 );

  // check that interpolation of polynomial is accurate to machine precision
  for( size_t i=0; i<n-1; ++i )
  {
    if( std::abs(oldErrors2[ i ]) > 1e-11 )
      DUNE_THROW(Dune::GridError,"interpolation of polynomial of order " << polOrder << " not exact!");
  }

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

void runTest( const int refCount, const int steps, std::istream& gridfile )
{
  // create grid
  Dune::GridPtr< GridType > grid( std::to_string( GridType::dimension ) + "dgrid.dgf" );
  // Dune::GridPtr< GridType > grid( gridfile );
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
  Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
        discontinuousGalerkinSpace( gridPart );
  Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
        lagrangeDGSpaceA( gridPart );
#if HAVE_DUNE_LOCALFUNCTIONS
  Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpaceType, GridPartType, polOrder >
        lagrangeDGSpaceB( gridPart );
  Dune::Fem::DGLagrangeSpace< FunctionSpaceType, GridPartType,
      Dune::GaussLobattoPointSet
  > lobattoDGSpace( gridPart, polOrder );
  Dune::Fem::FixedOrderDGLagrangeSpace< FunctionSpaceType, GridPartType, polOrder,
      Dune::GaussLegendrePointSet
  > gaussDGSpace( gridPart );
#endif
  Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder >
        lagrangeSpaceA( gridPart );
#if HAVE_DUNE_LOCALFUNCTIONS
  Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType >
        lagrangeSpaceB( gridPart, polOrder );
  Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType,
      Dune::GaussLobattoPointSet
  > lobattoSpace( gridPart, polOrder );
#endif
  Dune::Fem::HierarchicLegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, polOrder >
        hLegendreSpace( gridPart );
  Dune::Fem::FourierDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >
        fourierSpace( gridPart, polOrder+1 );

#if 0
  // Test:
  {
    auto quadrature = lobattoDGSpace.quadrature(Dune::GeometryTypes::cube(3));
    double weight = 0;
    // std::cout << "---- Lobatto ------\n";
    for (auto &p : quadrature)
    {
      // std::cout << p.point() << " , " << p.weight() << " , " << p.localKey() << std::endl;
      weight += p.weight();
    }
    assert(std::abs(weight-1.)<1e-15);
  }
  {
    auto quadrature = gaussDGSpace.quadrature(Dune::GeometryTypes::cube(3));
    double weight = 0;
    // std::cout << "---- Gauss ------\n";
    for (auto &p : quadrature)
    {
      // std::cout << p.point() << " , " << p.weight() << " , " << p.localKey() << std::endl;
      weight += p.weight();
    }
    assert(std::abs(weight-1.)<1e-15);
    // std::cout << "--------------------\n";
  }
#endif
  // perform eoc loop
  eocLoop( *grid, steps,
           discontinuousGalerkinSpace,
#if HAVE_DUNE_LOCALFUNCTIONS
           lagrangeDGSpaceB,
           lobattoDGSpace,
           gaussDGSpace,
#endif
           lagrangeDGSpaceA,
           lagrangeSpaceA,
#if HAVE_DUNE_LOCALFUNCTIONS
           lagrangeSpaceB,
           lobattoSpace,
#endif
           hLegendreSpace,
           fourierSpace
         );
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

  std::ifstream cartesiangrid( std::to_string( GridType::dimension ) + "dgrid.dgf" );
  if( cartesiangrid.is_open() )
  {
    runTest( refCount, steps, cartesiangrid );
  }

  // same tests on a non-affine grid if the grid is capable of this
  if( Dune::Capabilities::isCartesian< GridType >::v ) return 0;

  const double coords[8][3] = {
            {  0.0,  0.0,  0.0 }, // p0
            {  1.5, -0.1,  0.1 }, // p1
            {  0.1,  0.5,  0.2 }, // p2
            {  1.3,  0.4, -0.1 }, // p3
            { -0.3,  0.1,  1.2 }, // p4
            {  1.5, -0.1,  0.9 }, // p5
            {  0.2,  0.7,  1.0 }, // p6
            {  0.7,  0.6,  0.6 }};// p7

  std::stringstream nonAffineGrid;
  nonAffineGrid << "DGF" << std::endl;
  nonAffineGrid << "Vertex" << std::endl;
  const int vx = GridType::dimension == 3 ? 8 : 4;
  for( int i=0; i<vx; ++i )
  {
    for( int d=0; d<GridType::dimension; ++d )
    {
      nonAffineGrid << coords[ i ][ d ] << " ";
    }
    nonAffineGrid << std::endl;
  }
  nonAffineGrid << "#" << std::endl;

  nonAffineGrid << "Cube" << std::endl;
  for( int i=0; i<vx; ++i )
    nonAffineGrid << i << " ";
  nonAffineGrid << std::endl << "#" << std::endl;
  nonAffineGrid << "Simplex" << std::endl;
  nonAffineGrid << "#" << std::endl;

  runTest( refCount, steps, nonAffineGrid );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e << std::endl;
  return 1;
}
