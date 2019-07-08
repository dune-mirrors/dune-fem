#include <config.h>

//#define WANT_CACHED_COMM_MANAGER 0
//#define DG_ONLY

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
#include <dune/fem/space/common/adaptationmanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/space/common/restrictprolongtuple.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>


using Dune::elements;
using Dune::Fem::gridFunctionAdapter;

static const int polOrder = 2;
static const int dgOrder = 2;

static const int dimension = 2;


typedef Dune::Fem::FunctionSpace< double, double, dimension, 1 > FunctionSpaceType;

#if HAVE_DUNE_ALUGRID
typedef Dune::ALUGrid< dimension, dimension, Dune::simplex, Dune::conforming > GridType;
//typedef Dune::ALUGrid< dimension, dimension, Dune::simplex, Dune::nonconforming > GridType;

//typedef Dune::Fem::AdaptiveLeafGridPart< GridType, Dune::InteriorBorder_Partition > GridPartType;
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;
//typedef Dune::Fem::DGAdaptiveLeafGridPart< GridType > DGGridPartType;

typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, dgOrder > DGSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< DGSpaceType > DGFunctionType;

#ifdef DG_ONLY
typedef DGSpaceType     LagrangeSpaceType;
typedef DGFunctionType  LagrangeFunctionType;
#else
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, polOrder > LagrangeSpaceType;
typedef Dune::Fem::AdaptiveDiscreteFunction< LagrangeSpaceType > LagrangeFunctionType;
#endif

typedef Dune::Fem::RestrictProlongDefaultTuple< LagrangeFunctionType, DGFunctionType > RestrictProlongType;
typedef Dune::Fem::AdaptationManager< GridType, RestrictProlongType > AdaptationManagerType;
static const int maxLevel = 3 * Dune::DGFGridInfo< GridType > :: refineStepsForHalf();
#else
static const int maxLevel = 3;
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
  auto y = element.geometry().center();
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
  std::stringstream gridin;
  Dune::FieldVector<int, dimension> lower(0), upper(1), cells(8);
  gridin << "DGF" << std::endl;
  gridin << "INTERVAL" << std::endl;
  gridin << lower << std::endl;
  gridin << upper << std::endl;
  gridin << cells << std::endl;
  gridin << "#" << std::endl;

  Dune::GridPtr< GridType > grid( gridin );
  // if grid is not load balanced then one core has no elements and the test
  // fails. Needs to be checked.
  grid->loadBalance();

  GridPartType   gridPart( *grid );
  LagrangeSpaceType lagrangeSpace( gridPart );
  DGSpaceType dgSpace( gridPart );

  LagrangeFunctionType lagrangeSolution( "lagrange-solution", lagrangeSpace );
  interpolate( gridFunctionAdapter( ExactSolution(), gridPart, polOrder ), lagrangeSolution );

  DGFunctionType dgSolution( "dg-solution", dgSpace );
  interpolate( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder ), dgSolution  );

  // compute initial error
  Dune::Fem::L2Norm< GridPartType > l2Norm( gridPart );

  // compute initial error
  const double initialLagrangeErrorBefore =
    l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution );
  const double initialDGErrorBefore =
    l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution );

  // redistribute the grid before initializing the
  RestrictProlongType restrictProlong( lagrangeSolution, dgSolution );
  AdaptationManagerType adaptManager( *grid, restrictProlong );

  // compute initial error
  const double initialLagrangeError =
    l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution, Dune::Partitions::all );
  const double initialDGError =
    l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution, Dune::Partitions::all );

  {
    const double lagrangeErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution, Dune::Partitions::ghost );
    const double dgErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution, Dune::Partitions::ghost );

    // substract error on ghost cells from initial error
    const double lagrangeError = std::sqrt( initialLagrangeError*initialLagrangeError - lagrangeErrorGhosts*lagrangeErrorGhosts );
    if( std::abs( initialLagrangeErrorBefore - lagrangeError ) > 1e-12 )
      DUNE_THROW(Dune::ParallelError,"Lagrange error does not match after load balance");
    const double dgError = std::sqrt( initialDGErrorBefore*initialDGErrorBefore - dgErrorGhosts*dgErrorGhosts );
    if( std::abs( initialDGErrorBefore - dgError ) > 1e-12 )
      DUNE_THROW(Dune::ParallelError,"DG error does not match after load balance");
  }


  if( grid->comm().rank() == 0 )
    std::cout << "initial L2 error: " << std::scientific << std::setprecision( 12 ) << initialLagrangeError << "    " << initialDGError << std::endl;

  const bool doOutput = argc > 1 ;

  Dune::Fem::VTKIO< GridPartType > output( gridPart, Dune::VTK::nonconforming );
  if( doOutput )
  {
    output.addVertexData( lagrangeSolution );
    output.addVertexData( dgSolution );
    output.write( "balladapt-initial" );
  }

  double time = 0;

  // perform initial refinement
  for( int i = 0; i < maxLevel; ++i )
  {
    for( const auto &element : elements( gridPart ) )
      grid->mark( marking( element, time ), element );

    // do adaptation and loadBalance
    adaptManager.adapt();
    //adaptManager.loadBalance();
  }

  // let the ball rotate
  for( int nr = 0; time < 2.0*M_PI; time += 0.1, ++nr )
  {
    for( const auto &element : elements( gridPart ) )
      grid->mark( marking( element, time ), element );

    // do adaptation and loadBalance
    adaptManager.adapt();
    //adaptManager.loadBalance();

    if( doOutput )
    {
      output.write( "balladapt-" + std::to_string( nr ) );
    }
  }

  // coarsen up to macro level again
  for( int i = 0; i < maxLevel; ++i )
  {
    if( grid->maxLevel() == 0 )
      break;

    for( const auto &element : elements( gridPart ) )
      grid->mark( -1, element );

    // do adaptation and loadBalance
    adaptManager.adapt();
    //adaptManager.loadBalance();
  }

  if( doOutput )
  {
    output.write( "balladapt-final" );
  }
  if( grid->maxLevel() > 0 )
    DUNE_THROW( Dune::GridError, "Unable to coarsen back to macro grid" );

  const double finalLagrangeError = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution, Dune::Partitions::all );
  const double finalDGError = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution, Dune::Partitions::all );

  if( grid->comm().rank() == 0 )
    std::cout << "final L2 error:   " << std::scientific << std::setprecision( 12 ) << finalLagrangeError << "    " << finalDGError << std::endl;

  if( grid->leafGridView().ghostSize( 0 ) > 0 )
  {
    const double lagrangeErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution, Dune::Partitions::ghost );
    const double dgErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution, Dune::Partitions::ghost );

    // compare errors with "correct" errors
    interpolate( gridFunctionAdapter( ExactSolution(), gridPart, polOrder ), lagrangeSolution );
    interpolate( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder ), dgSolution  );

    const double trueLagrangeErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, polOrder+1 ), lagrangeSolution, Dune::Partitions::ghost );
    const double trueDGErrorGhosts = l2Norm.distance( gridFunctionAdapter( ExactSolution(), gridPart, dgOrder+1 ), dgSolution, Dune::Partitions::ghost );

    if( grid->comm().rank() == 0 )
    {
      std::cout << "Ghosts (projected - observed) "
                << std::abs( lagrangeErrorGhosts - trueLagrangeErrorGhosts ) << " "
                << std::abs( dgErrorGhosts - trueDGErrorGhosts ) << std::endl;
    }

    if( std::abs( lagrangeErrorGhosts - trueLagrangeErrorGhosts ) > 1e-12 || lagrangeErrorGhosts < 1e-10 )
      DUNE_THROW(Dune::ParallelError,"Lagrange error on ghost cells does not match after load balance");
    if( std::abs( dgErrorGhosts - trueDGErrorGhosts ) > 1e-12 || dgErrorGhosts < 1e-10 )
      DUNE_THROW(Dune::ParallelError,"DG error on ghost cells does not match after load balance");
  }

  return ((std::abs( initialLagrangeError - finalLagrangeError) < 1e-10) && (std::abs( initialDGError - finalDGError ) < 1e-10) ? 0 : 1);
#else // #if HAVE_DUNE_ALUGRID
  return 0;
#endif // #else // #if HAVE_DUNE_ALUGRID
}
catch( const Dune::Exception &exception )
{
  std::cerr << exception << std::endl;
  return 1;
}
