#include <config.h>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

//#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/gridpart/leafgridpart.hh>

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


int main ( int argc, char *argv[] )
{
  Dune::Fem::MPIManager::initialize( argc, argv );

  // grid type
  typedef Dune::YaspGrid< 2 > GridType;

  // create grid
  std::istringstream dgf( dgfUnitCube( GridType::dimensionworld, 4 ) );
  Dune::GridPtr< GridType > gridPtr( dgf );
  GridType& grid = *gridPtr;

  typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
  GridPartType gridPart( grid );

  typedef Dune::Fem::FunctionSpace<double,double,GridType::dimension, 1 >  FunctionSpaceType;
  typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 2 >
    DiscreteSpaceType ;

  DiscreteSpaceType space( gridPart );
  std::shared_ptr< DiscreteSpaceType > weakPtr1 = Dune::Fem::referenceToSharedPtr( space );
  std::shared_ptr< DiscreteSpaceType > weakPtr2 = Dune::Fem::referenceToSharedPtr( space );

  if( weakPtr1.operator->() != weakPtr2.operator->() )
    DUNE_THROW(Dune::GridError,"referenceToSharedPtr not working correctly");

  if( weakPtr1.use_count() != 2 )
    DUNE_THROW(Dune::GridError,"referenceToSharedPtr not working correctly");

  std::shared_ptr< DiscreteSpaceType > spcPtr( new DiscreteSpaceType(gridPart) );
  // create a reference
  DiscreteSpaceType& spc = *spcPtr;

  // now if the reference is used to create another shared ptr it should point
  // to the original shared ptr reference
  std::shared_ptr< DiscreteSpaceType > spcPtr2 =
    Dune::Fem::referenceToSharedPtr( spc );

  if( spcPtr2.operator->() != spcPtr.operator->() )
    DUNE_THROW(Dune::GridError,"referenceToSharedPtr not working correctly");

  if( spcPtr.use_count() != 2 )
    DUNE_THROW(Dune::GridError,"referenceToSharedPtr not working correctly");

  return 0;
}
