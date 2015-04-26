#include <config.h>
#include <iostream>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/finitevolume.hh>

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction.hh>
#endif

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/combinedfunction.hh>

#include <dune/fem/misc/mpimanager.hh>

typedef Dune:: GridSelector::GridType HGridType;
typedef Dune::Fem::DGAdaptiveLeafGridPart< HGridType > GridPartType;
typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, HGridType::dimension+2 > FunctionSpaceType;
typedef Dune::Fem::FiniteVolumeSpace< FunctionSpaceType, GridPartType, 0 > DiscreteFunctionSpaceType;

template <class DiscreteFunction, class OtherDiscreteFunction>
void checkFunction( DiscreteFunction& df, const OtherDiscreteFunction& other )
{
  typedef typename DiscreteFunction :: DofType DofType;
  std::fill( df.dbegin(), df.dend(), DofType( 0 ) );

  df += other;
  df -= other;

  df.assign( other );

  df *= 2.0;
  df /= 2.0;
}

// main program
int main(int argc, char ** argv)
{
  Dune::Fem::MPIManager :: initialize( argc, argv );
  try
  {
    HGridType &grid = Dune::Fem::TestGrid :: grid();

    GridPartType gridPart( grid );
    // add check for grid width
    std::cout << "Grid width: "
      << Dune::Fem::GridWidth :: calcGridWidth( gridPart ) << std::endl;

    DiscreteFunctionSpaceType space( gridPart );

    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > ref ("ref", space);

    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > adf ("adaptive", space);
    checkFunction( adf, ref );

    Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > istldf ("adaptive", space);
    checkFunction( istldf, ref );

    checkFunction( adf, istldf );

    return 0;
  }
  catch( const Dune::Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
