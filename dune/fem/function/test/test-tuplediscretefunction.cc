#include <config.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/lagrange.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/tuplediscretefunction.hh>

#include <dune/fem/misc/mpimanager.hh>

typedef Dune::GridSelector::GridType HGridType;
typedef Dune::Fem::AdaptiveLeafGridPart< HGridType > GridPartType;

typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, 1 > BaseFunctionSpaceType;

typedef Dune::Fem::DiscontinuousGalerkinSpace< BaseFunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpace1;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< BaseFunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpace2;
typedef Dune::Fem::DiscontinuousGalerkinSpace< BaseFunctionSpaceType, GridPartType, 1 > DiscreteFunctionSpace3;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< BaseFunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpace4;

typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace1 > DiscreteFunction1;
typedef Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpace2 > DiscreteFunction2;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace3 > DiscreteFunction3;
typedef Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpace4 > DiscreteFunction4;

typedef Dune::Fem::TupleDiscreteFunction< DiscreteFunction1, DiscreteFunction2, DiscreteFunction3, DiscreteFunction4 >
DiscreteFunctionType;

typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

/*
template <class DiscreteFunction, class OtherDiscreteFunction>
void checkFunction( DiscreteFunction& df, OtherDiscreteFunction& other )
{
  std::cout << "Checking (" << df.name() << "," << other.name()
            << "), size = ("<<df.size()<< "," << other.size() << ")....";
  typedef typename DiscreteFunction :: DofType DofType;

  // fill df with zeros
  df.clear();
  other.clear();

  // fill df with zeros
  std::fill( df.dbegin(), df.dend(), DofType( 0 ) );

  df += other;
  df -= other;

  // fill df with zeros
  for( auto&& dof : dofs(df) )
    dof=static_cast<DofType>(0);

  // fill with increasing values
  int cont(0);
  for( const auto& entity : entities(df) )
  {
    auto lf = df.localFunction( entity );
    lf.clear();
    for( auto i=0; i<lf.numDofs(); ++i,++cont )
      lf[ i ] = static_cast<DofType>(cont);
  }

  // check block access
  const size_t localBlockSize = DiscreteFunctionSpaceType::localBlockSize;
  const size_t numBlocks      = df.blocks();
  if( df.size() / localBlockSize != numBlocks )
    DUNE_THROW(Dune::InvalidStateException,"number of blocks not correct!");

  auto dfDofIt(df.dbegin());
  for(size_t i=0;i<numBlocks;++i)
    for(size_t j=0;j<localBlockSize;++j,++dfDofIt)
      if( std::abs( df.dofVector()[i][j] - *dfDofIt ) > 1e-12 )
        DUNE_THROW(Dune::InvalidStateException,"Block access did not work");

  // copy to std::vector, sometimes needed for solver interfaces
  std::vector< DofType > vec( df.size() );
  std::copy( df.dbegin(), df.dend(), vec.begin() );

  // check copy constructor
  DiscreteFunction copydf( df );
  if( ! (copydf == df) )
  {
    assert( false );
    DUNE_THROW(Dune::InvalidStateException,"Copying did not work");
  }

  df.dofVector()[0] *= 1.0;
  df.dofVector()[0][0] = 1.0;
  df.dofVector()[0][HGridType::dimension-1] = 1.0;

  df.assign( other );

  df *= 2.0;
  df /= 2.0;

  df.axpy( 0.5, other );

  df.enableDofCompression();

  std::stringstream stream;
  Dune::Fem::StandardOutStream out( stream );
  df.write( out );
  Dune::Fem::StandardInStream in( stream );
  df.read( in );

  df.scalarProductDofs( other );
  df.scalarProductDofs( df );
  other.scalarProductDofs( df );

  std::cout << "done!" << std::endl;

  typedef Dune::Fem::RestrictProlongDefault<DiscreteFunction> RPDefaultType;
  RPDefaultType rp( df );
  rp.setFatherChildWeight(Dune::DGFGridInfo< HGridType >::refineWeight());

  //typedef Dune::Fem::AdaptationManager< HGridType, RPDefaultType > AdaptationManagerType;
  //AdaptationManagerType adop(grid,rp);
}
*/

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

    DiscreteFunctionType df( "ref", space );
    std::cout << "dofs = " << df.size() << std::endl;

    // refine grid
    Dune::Fem::GlobalRefine::apply( grid, 1 );

    std::cout << "dofs = " << df.size() << std::endl;

    return 0;
  }
  catch( const Dune::Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
