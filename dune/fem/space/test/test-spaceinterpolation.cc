#include <config.h>

#include <cmath>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/classname.hh>
#include <dune/common/ftraits.hh>
#include <dune/common/hybridutilities.hh>

#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/common/gridfunctionadapter.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/interpolate.hh>
#include <dune/fem/test/exactsolution.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/misc/h1norm.hh>

#include <dune/fem/space/brezzidouglasmarini.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/padaptivespace.hh>
#include <dune/fem/space/rannacherturek.hh>

#include <dune/fem/space/test/checklocalinterpolation.hh>


// dgfUnitCube
// -----------

inline static std::string dgfUnitCube ( int dimWorld, int cells )
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
  dgf += "\n#\n";
  return dgf;
}



// Real
// ----

template< class T >
using Real = typename Dune::FieldTraits< typename T::RangeFieldType >::real_type;



// ErrorTuple
// ----------

template< class T >
struct ErrorTuple;

template< class... T >
struct ErrorTuple< std::tuple< T... > >
{
  typedef std::tuple< std::pair< Real< T >, Real< T > >... > Type;
};



// Type Definitions
// ----------------

typedef Dune::GridSelector::GridType GridType;
typedef Dune::Fem::LeafGridPart< GridType > GridPartType;

static const int dimRange = GridPartType::dimensionworld;

typedef Dune::Fem::GridFunctionSpace< GridPartType, Dune::FieldVector< typename GridPartType::ctype, dimRange > > FunctionSpaceType;

typedef std::tuple<
  Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType >,
  Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 0 >,
  Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 1 >,
  Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 2 >,
  Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 1 >,
  Dune::Fem::LagrangeDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 2 >,
  Dune::Fem::LegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 1 >,
  Dune::Fem::LegendreDiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 2 >,
#if HAVE_DUNE_LOCALFUNCTIONS
  Dune::Fem::BrezziDouglasMariniSpace< FunctionSpaceType, GridPartType, 1 >,
  Dune::Fem::BrezziDouglasMariniSpace< FunctionSpaceType, GridPartType, 2 >,
  // Dune::Fem::LagrangeSpace< FunctionSpaceType, GridPartType >,
  // Dune::Fem::RannacherTurekSpace< FunctionSpaceType, GridPartType >,
#endif // #if HAVE_DUNE_LOCALFUNCTIONS
  Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 1 >,
  Dune::Fem::LagrangeDiscreteFunctionSpace< FunctionSpaceType, GridPartType, 2 >
  > DiscreteFunctionSpacesType;

typedef ErrorTuple< DiscreteFunctionSpacesType >::Type ErrorTupleType;



// algorithm
// ---------

template< class DiscreteFunctionSpace >
std::pair< Real< DiscreteFunctionSpace >, Real< DiscreteFunctionSpace > >
algorithm ( typename DiscreteFunctionSpace::GridPartType &gridPart )
{
  DiscreteFunctionSpace space( gridPart );
  Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace > u( "solution", space );

  // interpolate a function
  Dune::Fem::ExactSolution< typename DiscreteFunctionSpace::FunctionSpaceType > uExact;
  const auto uGridExact = gridFunctionAdapter( "exact solution", uExact, gridPart, 3 );
  interpolate( uGridExact, u );

  checkLocalInterpolation( space );

  Dune::Fem::L2Norm< GridPartType > l2norm( gridPart );
  Dune::Fem::H1Norm< GridPartType > h1norm( gridPart );

  return std::make_pair( l2norm.distance( uGridExact, u ), h1norm.distance( uGridExact, u ) );
}



// eoc
// ---

template< class T, class U >
std::pair< T, U > eoc ( const std::pair< T, U > &e_old, const std::pair< T, U > &e_new )
{
  using std::log;
  return std::make_pair( log( e_old.first / e_new.first ) / log( 2 ), log( e_old.second / e_new.second ) / log( 2 ) );
}


// print
// -----

template< class L2Error, class L2Eoc, class H1Error, class H1Eoc >
void print ( L2Error &&l2Error, L2Eoc &&l2Eoc, H1Error &&h1Error, H1Eoc &&h1Eoc )
{
  std::cout << std::setw( 12 ) << std::setprecision( 6 ) << l2Error;
  std::cout << std::setw( 12 ) << std::setprecision( 2 ) << l2Eoc;
  std::cout << std::setw( 12 ) << std::setprecision( 6 ) << h1Error;
  std::cout << std::setw( 12 ) << std::setprecision( 2 ) << h1Eoc;
  std::cout << std::endl;
}



// main
// ----

int main ( int argc, char **argv )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // construct unit cube
  typedef typename Dune::GridSelector::GridType GridType;
  std::istringstream dgf( dgfUnitCube( GridType::dimensionworld, 4 ) );
  Dune::GridPtr< GridType > grid( dgf );

  // create leaf grid part
  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( *grid );

  auto indices = std::make_index_sequence< std::tuple_size< DiscreteFunctionSpacesType >::value >();

  std::array< ErrorTupleType, 4 > errors;
  for( ErrorTupleType &e : errors )
  {
    Dune::Hybrid::forEach( indices, [ &gridPart, &e ] ( auto &&idx ) {
        const std::size_t i = std::decay_t< decltype( idx ) >::value;
        std::get< i >( e ) = algorithm< std::tuple_element_t< i, DiscreteFunctionSpacesType > >( gridPart );
      } );
    grid->globalRefine(1);
  }

  Dune::Hybrid::forEach( indices, [ &errors ] ( auto &&idx ) {
      const std::size_t i = std::decay_t< decltype( idx ) >::value;
      std::cout << ">>> Testing " << Dune::className< std::tuple_element_t< i, DiscreteFunctionSpacesType > >() << ":" << std::endl;
      std::cout << std::endl;
      print( "L2 Error", "L2 EOC", "H1 Error", "H1 EOC" );
      print( std::get< i >( errors[ 0 ] ).first, "---", std::get< i >( errors[ 0 ] ).second, "---" );
      for( std::size_t j = 1; j < errors.size(); ++j )
      {
        auto eocs = eoc( std::get< i >( errors[ j-1 ] ), std::get< i >( errors[ j ] ) );
        print( std::get< i >( errors[ j ] ).first, eocs.first, std::get< i >( errors[ j ] ).second, eocs.second );
      }
      std::cout << std::endl;
    } );

  return 0;
}