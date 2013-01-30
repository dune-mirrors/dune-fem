#include <config.h>

#include <iostream>
#include <sstream>
#include <string>

#include <dune/common/exceptions.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/test/testgrid.hh>

#include "../space.hh"


static const int order = ORDER;

typedef Dune::GridSelector::GridType GridType;
typedef Dune::Fem::AdaptiveLeafGridPart< GridType > GridPartType;

static const int dimDomain = GridPartType::dimensionworld;
static const int dimRange = 1;

typedef Dune::Fem::FunctionSpace< GridPartType::ctype, double, dimDomain, dimRange > FunctionSpaceType;
typedef Dune::Fem::FourierDiscreteFunctionSpace< FunctionSpaceType, GridPartType, order > DiscreteFunctionSpaceType;



template< class GridType >
void test ( GridType &grid )
{
  GridPartType gridPart( grid );
  DiscreteFunctionSpaceType discreteFunctionSpace( gridPart );

  typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > DiscreteFunctionType;
  DiscreteFunctionType solution( "solution", discreteFunctionSpace );
  solution.clear();

  typedef typename GridPartType::template Codim< 0 >::IteratorType IteratorType;
  typedef typename IteratorType::Entity EntityType;

  const IteratorType end = gridPart.template end< 0 >();
  for( IteratorType it = gridPart.template begin< 0 >(); it != end; ++it )
  {
    const EntityType &entity = *it;

    typedef typename DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
    BasisFunctionSetType basisFunctionSet = discreteFunctionSpace.basisFunctionSet( entity );
    if( !basisFunctionSet.size() == discreteFunctionSpace.size() )
      DUNE_THROW( Dune::InvalidStateException, "Sizes do not match." );
  }
}



int main ( int argc, char **argv )
try
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  Dune::Fem::Parameter::append( argc, argv );
  Dune::Fem::Parameter::append( "parameter" );

  const int refCount = Dune::Fem::Parameter::getValue< int >( "startLevel", 0 );

  GridType &grid = Dune::Fem::TestGrid::grid();
  const int refineStepsForHalf = Dune::Fem::TestGrid::refineStepsForHalf();
  for( int count = 0; count < refCount; ++count )
    Dune::Fem::GlobalRefine::apply( grid, refineStepsForHalf );

  test( grid );

  return 0;
}
catch( const Dune::Exception &e )
{
  std::cerr << e.what() << std::endl;
  return 0;
}
catch( const std::exception &e )
{
  std::cerr << e.what() << std::endl;
  return 1;
}
