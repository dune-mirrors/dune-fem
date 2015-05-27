#include <config.h>
#include <iostream>

#include <dune/fem/test/testgrid.hh>
#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/space/finitevolume.hh>
#include <dune/fem/space/discontinuousgalerkin.hh>
#include <dune/fem/space/common/adaptmanager.hh>

#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/storage/vector.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/function/vectorfunction/managedvectorfunction.hh>
#include <dune/fem/function/blockvectordiscretefunction.hh>
#include <dune/fem/function/blockvectors/referenceblockvector.hh>
#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction.hh>
#endif

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/combinedfunction.hh>

#include <dune/fem/space/common/adaptmanager.hh>
#include <dune/fem/misc/mpimanager.hh>

typedef Dune:: GridSelector::GridType HGridType;
typedef Dune::Fem::DGAdaptiveLeafGridPart< HGridType > GridPartType;
typedef Dune::Fem::FunctionSpace< double, double, HGridType::dimension, HGridType::dimension+2 > FunctionSpaceType;
typedef Dune::Fem::DiscontinuousGalerkinSpace< FunctionSpaceType, GridPartType, 2 > DiscreteFunctionSpaceType;

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

  // fill with increasing values
  int cont(0);
  for( const auto entity : df.space() )
  {
    auto lf = df.localFunction( entity );
    lf.clear();
    for( int i=0; i<lf.numDofs(); ++i,++cont )
      lf[ i ] = static_cast<DofType>(cont);
  }

  // check block access
  const std::size_t localBlockSize(DiscreteFunctionSpaceType::localBlockSize);
  const std::size_t numBlocks(df.size()/localBlockSize);
  typename DiscreteFunction::DofIteratorType dfDofIt(df.dbegin());
  for(std::size_t i=0;i!=numBlocks;++i)
    for(std::size_t j=0;j!=localBlockSize;++j,++dfDofIt)
      if((*df.block(i))[j]!=*dfDofIt)
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

  (*df.block( 0 )) *= 1.0;
  (*df.block( 0 ))[ 0 ] = 1.0;
  (*df.block( 0 ))[ HGridType::dimension-1 ] = 1.0;

  df.assign( other );

  df *= 2.0;
  df /= 2.0;

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
    std::cout << "dofs = " << adf.size() << std::endl;

    checkFunction( adf, ref );

    std::vector< double > advec( space.size() );
    Dune::Fem::AdaptiveDiscreteFunction< DiscreteFunctionSpaceType > adfp ("pointer", space, advec.data() );
    checkFunction( adfp, ref );

    Dune::Fem::DynamicVector< double > vec( space.size() );
    typedef Dune::Fem::VectorDiscreteFunction< DiscreteFunctionSpaceType, Dune::Fem :: DynamicVector< double > > VectorDiscreteFunctionType;
    VectorDiscreteFunctionType vdf ("vector", space, vec);
    checkFunction( vdf, ref );

    std::vector<double> vec1( space.size() );
    typedef Dune::Fem::VectorDiscreteFunction< DiscreteFunctionSpaceType, std::vector<double> > StdVectorDiscreteFunctionType;
    StdVectorDiscreteFunctionType vdf1 ("vector", space, vec1);
    checkFunction( vdf1, ref );

    Dune::Fem::ManagedDiscreteFunction< VectorDiscreteFunctionType > mdf ("managed", space);
    checkFunction( mdf, ref );

    typedef Dune::Fem::ReferenceBlockVector< FunctionSpaceType::RangeFieldType, DiscreteFunctionSpaceType::localBlockSize >
      BlockVectorType;
    Dune::Fem::BlockVectorDiscreteFunction< DiscreteFunctionSpaceType, BlockVectorType > bdf( "block", space );
    checkFunction( bdf, ref );

#if HAVE_DUNE_ISTL
    Dune::Fem::ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceType > istldf ("istl", space);
    checkFunction( istldf, ref );

    checkFunction( adf, istldf );
    checkFunction( vdf, istldf );
#endif

#if HAVE_PETSC
    Dune::Fem::PetscDiscreteFunction< DiscreteFunctionSpaceType > petscdf ("petsc", space);
    checkFunction( petscdf, ref );
    checkFunction( petscdf, vdf );
#endif

    // refine grid
    Dune::Fem::GlobalRefine::apply( grid, 1 );

    std::cout << "dofs = " << adf.size() << std::endl;
    std::cout << "dofs = " << istldf.size() << std::endl;
    checkFunction( adf, istldf );
    checkFunction( bdf, istldf );
    checkFunction( istldf, ref );
    checkFunction( mdf, ref );

    return 0;
  }
  catch( const Dune::Exception& e )
  {
    std :: cerr << e.what() << std :: endl;
    return 1;
  }
}
