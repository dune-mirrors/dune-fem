#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>

#include <dune/common/dynvector.hh>

#include <dune/fem/gridpart/idgridpart.hh>
#include <dune/fem/gridpart/leafgridpart.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/space/lagrange.hh>

#include <dune/fem/operator/linear/test/checklinearoperator.hh>

typedef Dune::GridSelector::GridType GridType;

typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > SpaceType;

typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
typedef Dune::Fem::IdGridPart< GridPartType > TestGridPartType;

typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, 1 > DiscreteSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, TestGridPartType, 1 > DiscreteTestSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, 2 > P2DiscreteSpaceType;

#if HAVE_DUNE_ISTL
// the UMFPack implementation in dune-istl is inflexible and only works with
// FieldMartrix. Here, we test MyBlock and therefore disable the directsolver
// for AMG.
#define DISABLE_AMG_DIRECTSOLVER 1
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>

// emulate block type being derived from Dune::FieldVector (e.g. in OPM)
namespace Dune {
  template <class K, int dim>
  struct MyBlock : public Dune::FieldVector< K, dim >
  {
    typedef Dune::FieldVector< K, dim > BaseType;
    using BaseType :: operator =;
  };

  template <class K, int dim>
  struct FieldTraits< MyBlock<  K, dim > > : public FieldTraits< typename MyBlock<  K, dim > :: BaseType >
  {};
}

template< class DSpace, class RSpace >
struct ISTLOperator
{
  typedef Dune::Fem::ISTLLinearOperator<
    Dune::Fem::ISTLBlockVectorDiscreteFunction< DSpace, Dune::MyBlock< typename DSpace::RangeFieldType, DSpace::localBlockSize > >,
    Dune::Fem::ISTLBlockVectorDiscreteFunction< RSpace, Dune::MyBlock< typename RSpace::RangeFieldType, RSpace::localBlockSize > >
    > type;
};
#endif

#if HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>

template< class DSpace, class RSpace >
struct PetscOperator
{
  typedef Dune::Fem::PetscLinearOperator<
    Dune::Fem::PetscDiscreteFunction< DSpace >,
    Dune::Fem::PetscDiscreteFunction< RSpace >
    > type;
};
#endif

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>

template< class DSpace, class RSpace >
struct SparseRowOperator
{
  typedef Dune::Fem::SparseRowLinearOperator<
    Dune::Fem::AdaptiveDiscreteFunction< DSpace >,
    Dune::Fem::AdaptiveDiscreteFunction< RSpace >
    > type;
};


template <template <class,class> class LinearOperator>
void algorithm( GridType& grid, const bool petsc = false )
{
  GridPartType gridPart( grid );
  TestGridPartType testGridPart( gridPart );

  DiscreteSpaceType space( gridPart );
  DiscreteTestSpaceType testSpace( testGridPart );
  P2DiscreteSpaceType p2Space( gridPart );

  {
    // check for same space
    typename LinearOperator< DiscreteSpaceType, DiscreteSpaceType >::type
    linOp( "linOp", space, space );

    std::vector< std::pair< int, int > > permutation =
    {
      std::make_pair( 0, 0 ),
      std::make_pair( 1, 1 ),
      std::make_pair( 2, 2 ),
      std::make_pair( 3, 3 )
    };
    checkLinearOperator( linOp, diagonalRange( space, space ), permutation );
  }

  {
    // check for same space size but different grids
    typename LinearOperator< DiscreteSpaceType, DiscreteTestSpaceType >::type
    linOp( "linOp", space, testSpace );

    std::vector< std::pair< int, int > > permutation =
    {
      std::make_pair( 0, 0 ),
      std::make_pair( 1, 1 ),
      std::make_pair( 2, 2 ),
      std::make_pair( 3, 3 )
    };
    checkLinearOperator( linOp, diagonalRange( space, testSpace ), permutation );
  }

  // the following does not work with PETSc
  if( petsc ) return ;

  {
    // check for different space sizes, but same grid
    typename LinearOperator< DiscreteSpaceType, P2DiscreteSpaceType >::type
    linOp( "linOp", space, p2Space );

    std::vector< std::pair< int, int > > permutation =
    {
      std::make_pair( 0, 0 ),
      std::make_pair( 2, 1 ),
      std::make_pair( 6, 2 ),
      std::make_pair( 8, 3 )
    };
    checkLinearOperator( linOp, diagonalRange( space, p2Space ), permutation );
  }

  {
    // check for different space sizes, but same grid
    typename LinearOperator< DiscreteTestSpaceType, P2DiscreteSpaceType >::type
    linOp( "linOp", testSpace, p2Space );

    std::vector< std::pair< int, int > > permutation =
    {
      std::make_pair( 0, 0 ),
      std::make_pair( 2, 1 ),
      std::make_pair( 6, 2 ),
      std::make_pair( 8, 3 )
    };
    checkLinearOperator( linOp, diagonalRange( testSpace, p2Space ), permutation );
  }
}


// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  // initialize MPI manager and PETSc
  Dune::Fem::MPIManager::initialize( argc, argv );

  // GridType grid( {1, 1}, {{2, 2}} );
  std::stringstream gridfile;
  gridfile << "DGF" << std::endl;
  gridfile << "Interval" << std::endl;
  Dune::FieldVector< double, GridType::dimension> lower( 0 );
  Dune::FieldVector< double, GridType::dimension> upper( 1 );
  Dune::FieldVector< int,    GridType::dimension> length( 2 );

  gridfile << lower  << std::endl;
  gridfile << upper  << std::endl;
  gridfile << length << std::endl;
  gridfile << "#" << std::endl;

  Dune::GridPtr< GridType > gridPtr( gridfile );
  GridType& grid = *gridPtr;


  // SparseRowLinearOperator
  {
    algorithm< SparseRowOperator >( grid );
  }

#if HAVE_DUNE_ISTL
  // ISTLLinearOperator
  {
    algorithm< ISTLOperator >( grid );
  }
#endif

#if HAVE_PETSC
  // PetscOperator
  {
    algorithm< PetscOperator >( grid );
  }
#endif

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
