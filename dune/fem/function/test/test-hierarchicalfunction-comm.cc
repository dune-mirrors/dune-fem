#include <config.h>

#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/geometry/dimension.hh>

#include <dune/grid/common/partitionset.hh>
#include <dune/grid/common/rangegenerators.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/function/hierarchical.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/function/localfunction/const.hh>
#include <dune/fem/function/localfunction/temporary.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#include <dune/fem/space/discontinuousgalerkin/space.hh>
#include <dune/fem/space/discontinuousgalerkin/tuple.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/test/exactsolution.hh>


// dgfUnitCube
// -----------

inline static std::string dgfUnitCube ( int dimWorld, int cells, int overlap = 1 )
{
  std::string dgf = "DGF\nINTERVAL\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += " 0";
  dgf += "\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += " 1";
  dgf += "\n";
  for( int i = 0; i < dimWorld; ++i )
    dgf += (" " + std::to_string( cells ));
  dgf += "\n#\n\n";
  dgf += "GRIDPARAMETER\n";
  dgf += "OVERLAP " + std::to_string( overlap ) + "\n";
  dgf += "#\n";
  return dgf;
}



// FunctionSpace
// -------------

template< class GridPart, int dimRange >
using FunctionSpace = Dune::Fem::GridFunctionSpace< GridPart, Dune::Dim< dimRange > >;



// TaylorHoodDGSpace
// -----------------

template< class GridPart >
using VelocityDGSpace = Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpace< GridPart, GridPart::dimensionworld >, GridPart, 2 >;

template< class GridPart >
using PressureDGSpace = Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpace< GridPart, 1 >, GridPart, 1 >;

template< class GridPart >
using TaylorHoodDGSpace = Dune::Fem::TupleDiscontinuousGalerkinSpace< VelocityDGSpace< GridPart >, PressureDGSpace< GridPart > >;



// LagrangeSpace
// -------------

template< class GridPart >
using LagrangeSpace = Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpace< GridPart, GridPart::dimensionworld+1 >, GridPart, 2 >;



// equals
// ------

template< class T >
bool equals ( const T &u, const T &v )
{
  typedef std::decay_t< decltype( std::declval< const T & >().two_norm2() ) > real_type;

  T w( u );
  w -= v;
  return (w.two_norm2() < 128 * std::numeric_limits< real_type >::epsilon());
}



// interpolateOnly
// ---------------

template< class GridFunction, class DiscreteFunction >
void interpolateOnly ( const GridFunction &u, DiscreteFunction &v )
{
  v.clear();

  Dune::Fem::ConstLocalFunction< GridFunction > uLocal( u );
  Dune::Fem::TemporaryLocalFunction< typename DiscreteFunction::DiscreteFunctionSpaceType > vLocal( v.space() );
  Dune::Fem::LocalInterpolation< typename DiscreteFunction::DiscreteFunctionSpaceType > interpolation( v.space() );

  // iterate over selected partition
  for( const auto& entity : elements( v.gridPart(), Dune::Partitions::all ) )
  {
    // initialize u, v to entity
    auto uGuard = bindGuard( uLocal, entity );
    auto vGuard = bindGuard( vLocal, entity );

    // bind interpolation to entity
    auto iGuard = bindGuard( interpolation, entity );

    // perform local interpolation
    interpolation( uLocal, vLocal.localDofVector() );

    // write interpolation into global DoF vector
    v.setLocalDofs( entity, vLocal.localDofVector() );
  }
}



// performTest
// -----------

template< class DiscreteFunctionSpace >
void performTest ( const DiscreteFunctionSpace &dfSpace )
{
  // some stupid type defs
  typedef Dune::Fem::HierarchicalDiscreteFunction< DiscreteFunctionSpace > DiscreteFunctionType;

  typedef typename DiscreteFunctionSpace::GridPartType GridPartType;
  typedef typename DiscreteFunctionSpace::LocalBlockIndices BlockIndices;

  typedef typename DiscreteFunctionType::DofType DofType;

  // interpolate a function
  Dune::Fem::ExactSolution< FunctionSpace< GridPartType, GridPartType::dimensionworld+1 > > uExact;
  const auto uGridExact = gridFunctionAdapter( "exact solution", uExact, dfSpace.gridPart(), 3 );
  DiscreteFunctionType u( "solution", dfSpace );
  interpolateOnly( uGridExact, u );

  // copy discrete function and clear auxiliary dofs
  DiscreteFunctionType w( u );
  for( const auto &auxiliaryDof : dfSpace.auxiliaryDofs() )
    Dune::Hybrid::forEach( BlockIndices(), [ &w, &auxiliaryDof ] ( auto &&j ) { w.dofVector()[ auxiliaryDof ][ j ] = DofType( 0 ); } );

  // make sure u and w differ
  if( equals( u.dofVector().array(), w.dofVector().array() ) )
    DUNE_THROW( Dune::InvalidStateException, "Unique representation does not differ from consistent representation." );

  // communicate w
  w.communicate();

  // now u and w should not differ
  if( !equals( u.dofVector().array(), w.dofVector().array() ) )
    DUNE_THROW( Dune::InvalidStateException, "Functions differ after communication." );
}



// main
// ----

int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // construct unit cube
  typedef typename Dune::GridSelector::GridType GridType;
  std::istringstream dgf( dgfUnitCube( GridType::dimensionworld, 16 ) );
  Dune::GridPtr< GridType > grid( dgf );

  // create leaf grid part
  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  // test Taylor-Hood DG space
  TaylorHoodDGSpace< GridPartType > taylorHoodDGSpace( gridPart );
  performTest( taylorHoodDGSpace );

  // test Lagrange space
  LagrangeSpace< GridPartType > lagrangeSpace( gridPart );
  performTest( lagrangeSpace );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cout << "Exception: " << e << std::endl;
  return 1;
}
