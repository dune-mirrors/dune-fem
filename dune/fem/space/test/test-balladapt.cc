#include <config.h>

#include <cmath>

#include <iostream>

#include <dune/grid/common/rangegenerators.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif // #if HAVE_DUNE_ALUGRID

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/file/vtkio.hh>
#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/lagrange.hh>


using Dune::elements;
using Dune::Fem::gridFunctionAdapter;

static const int maxLevel = 8;
static const int polOrder = 2;

const char dgf[]
  = "DGF\n"
    "\n"
    "INTERVAL\n"
    "0  0\n"
    "1  1\n"
    "8  8\n"
    "#\n";


typedef Dune::Fem::FunctionSpace< double, double, 2, 1 > FunctionSpaceType;

#if HAVE_DUNE_ALUGRID
typedef Dune::ALUGrid< 2, 2, Dune::simplex, Dune::conforming > GridType;

typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > DiscreteFunctionSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;

typedef Dune::Fem::RestrictProlongDefault< DiscreteFunctionType > RestrictProlongType;
typedef Dune::Fem::AdaptationManager< GridType, RestrictProlongType > AdaptationManagerType;
#endif // #if HAVE_DUNE_ALUGRID



// ExactSolution
// -------------

struct ExactSolution
  : public Dune::Fem::Function< FunctionSpaceType, ExactSolution >
{
  typedef FunctionSpaceType::RangeType RangeType;
  typedef FunctionSpaceType::RangeFieldType RangeFieldType;
  typedef FunctionSpaceType::DomainType DomainType;

  void evaluate ( const DomainType &x, RangeType &ret ) const
  {
    ret[ 0 ] = std::sin( M_PI*x[ 0 ] )*std::cos( M_PI*x[ 1 ] );
  }
};




// marking
// -------

template< class Element >
int marking ( const Element &element, double time )
{
  Dune::FieldVector< double, 2 > y = element.geometry().center();
  y[ 0 ] -= 0.5 + 0.2*std::cos( time );
  y[ 1 ] -= 0.5 + 0.2*std::sin( time );
  if( y.two_norm2() < 0.1*0.1 )
    return (element.level() < maxLevel ? 1 : 0);
  else
    return -1;
}



// main
// ----

int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

#if HAVE_DUNE_ALUGRID
  std::istringstream gridin( dgf );
  Dune::GridPtr< GridType > grid( gridin );
  // bug: this call to load balance is necessary while it should not
  grid->loadBalance();

  GridPartType gridPart( *grid );
  DiscreteFunctionSpaceType dfSpace( gridPart );

  DiscreteFunctionType solution( "solution", dfSpace );
  interpolate( gridFunctionAdapter( ExactSolution(), gridPart ), solution );

  // compute initial error
  Dune::Fem::L2Norm< GridPartType > l2Norm( gridPart );
  const double initialError = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder ), solution );
  if( grid->comm().rank() == 0 )
    std::cout << "initial L2 error: " << std::scientific << std::setprecision( 12 ) << initialError << std::endl;

  Dune::Fem::VTKIO< GridPartType > output( gridPart );
  output.addVertexData( solution );
  output.write( "balladapt-initial" );

  RestrictProlongType restrictProlong( solution );
  AdaptationManagerType adaptManager( *grid, restrictProlong );
  adaptManager.loadBalance();

  double time = 0;

  // perform initial refinement
  for( int i = 0; i < maxLevel; ++i )
  {
    for( const auto &element : elements( gridPart ) )
      grid->mark( marking( element, time ), element );
    adaptManager.adapt();
    adaptManager.loadBalance();
  }

  // let the ball rotate
  for( int nr = 0; time < 2.0*M_PI; time += 0.1, ++nr )
  {
    for( const auto &element : elements( gridPart ) )
      grid->mark( marking( element, time ), element );
    adaptManager.adapt();
    adaptManager.loadBalance();
    output.write( "balladapt-" + std::to_string( nr ) );
  }

  // coarsen up to macro level again
  for( int i = 0; i < maxLevel; ++i )
  {
    if( grid->maxLevel() == 0 )
      break;

    for( const auto &element : elements( gridPart ) )
      grid->mark( -1, element );
    adaptManager.adapt();
    adaptManager.loadBalance();
  }
  output.write( "balladapt-final" );
  if( grid->maxLevel() > 0 )
    DUNE_THROW( Dune::GridError, "Unable to coarsen back to macro grid" );

  const double finalError = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder ), solution );
  if( grid->comm().rank() == 0 )
    std::cout << "final L2 error:   " << std::scientific << std::setprecision( 12 ) << finalError << std::endl;

  return (initialError - finalError < 1e-10 ? 0 : 1);
#else // #if HAVE_DUNE_ALUGRID
  return 77;
#endif // #else // #if HAVE_DUNE_ALUGRID
}
catch( const Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
