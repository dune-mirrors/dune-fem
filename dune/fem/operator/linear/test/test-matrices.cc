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

typedef Dune::GridSelector::GridType GridType;

typedef Dune::Fem::FunctionSpace< double, double, GridType::dimensionworld, 1 > SpaceType;

typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
typedef Dune::Fem::IdGridPart< GridPartType > TestGridPartType;

typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, 1 > DiscreteSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, TestGridPartType, 1 > DiscreteTestSpaceType;
typedef Dune::Fem::LagrangeDiscreteFunctionSpace< SpaceType, GridPartType, 2 > P2DiscreteSpaceType;



#if USE_ISTL && HAVE_DUNE_ISTL
#include <dune/fem/function/blockvectorfunction.hh>
#include <dune/fem/operator/linear/istloperator.hh>

template< class DSpace, class RSpace >
struct LinearOperator
{
  typedef Dune::Fem::ISTLLinearOperator<
    Dune::Fem::ISTLBlockVectorDiscreteFunction< DSpace >,
    Dune::Fem::ISTLBlockVectorDiscreteFunction< RSpace >
    > type;
};

#elif USE_PETSC && HAVE_PETSC
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/operator/linear/petscoperator.hh>

template< class DSpace, class RSpace >
struct LinearOperator
{
  typedef Dune::Fem::PetscLinearOperator<
    Dune::Fem::PetscDiscreteFunction< DSpace >,
    Dune::Fem::PetscDiscreteFunction< RSpace >
    > type;
};

#elif USE_EIGEN && HAVE_EIGEN
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/operator/linear/eigenoperator.hh>
#include <dune/fem/storage/eigenvector.hh>

template< class DSpace, class RSpace >
struct LinearOperator
{
  typedef Dune::Fem::EigenVector< double > DofVectorType;
  typedef Dune::Fem::EigenLinearOperator<
    Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< DSpace, DofVectorType > >,
    Dune::Fem::ManagedDiscreteFunction< Dune::Fem::VectorDiscreteFunction< RSpace, DofVectorType > >
    > type;
};

#else
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/linear/spoperator.hh>

template< class DSpace, class RSpace >
struct LinearOperator
{
  typedef Dune::Fem::SparseRowLinearOperator<
    Dune::Fem::AdaptiveDiscreteFunction< DSpace >,
    Dune::Fem::AdaptiveDiscreteFunction< RSpace >
    > type;
};

#endif


// DiagonalRange
// -------------

template< class DomainSpace, class RangeSpace >
struct DiagonalRange
{
  struct Iterator
  {
    typedef typename DomainSpace::IteratorType Iterator1;
    typedef typename DomainSpace::EntityType Entity1;
    typedef typename RangeSpace::IteratorType Iterator2;
    typedef typename RangeSpace::EntityType Entity2;
    Iterator ( std::pair< Iterator1, Iterator2 > it ) : it_( it ) {}

    void operator++ () { it_.first++; it_.second++; }

    std::pair< Entity1, Entity2 >
    operator* ()
    {
      return std::make_pair( *it_.first, *it_.second );
    }

    bool operator!= ( const Iterator &other ) const { return it_ != other.it_; }
    bool operator== ( const Iterator &other ) const { return it_ == other.it_; }

  private:
    std::pair< Iterator1, Iterator2 > it_;
  };

  typedef Iterator Iterator;

  DiagonalRange ( const DomainSpace &dSpace, const RangeSpace &rSpace )
    : dSpace_( dSpace ),
      rSpace_( rSpace )
  {}

  Iterator begin () const { return Iterator( std::make_pair( dSpace_.begin(), rSpace_.begin() ) ); }
  Iterator end () const { return Iterator( std::make_pair( dSpace_.end(), rSpace_.end() ) ); }

protected:
  const DomainSpace &dSpace_;
  const RangeSpace &rSpace_;
};


// diagonalRange
// -------------

template< class Space1, class Space2 >
DiagonalRange< Space1, Space2 > diagonalRange ( const Space1 &space1, const Space2 &space2 )
{
  return DiagonalRange< Space1, Space2 >( space1, space2 );
}


// RangeStencil
// ------------

template< class DomainSpace, class RangeSpace >
struct RangeStencil
  : public Dune::Fem::Stencil< DomainSpace, RangeSpace >
{
  typedef RangeStencil< DomainSpace, RangeSpace > ThisType;
  typedef Dune::Fem::Stencil< DomainSpace, RangeSpace > BaseType;

  template< class Range >
  RangeStencil ( const DomainSpace &dSpace, const RangeSpace &rSpace, const Range &range )
    : BaseType( dSpace, rSpace )
  {
    for( const auto &entry : range )
      BaseType::fill( entry.first, entry.second );
  }
};


// checkLinearOperator
// -------------------

template< class LinearOperator, class Range >
void checkLinearOperator ( LinearOperator &linOp, const Range &range, const std::vector< std::pair< int, int > > &permutation )
{
  // get type of domain and range function
  typedef typename LinearOperator::DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
  typedef typename LinearOperator::RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

  // type of localMatrix
  typedef typename LinearOperator::LocalMatrixType LocalMatrixType;

  // get spaces
  const DomainSpaceType &dSpace = linOp.domainSpace();
  const RangeSpaceType &rSpace = linOp.rangeSpace();

  // construct some domain and range functions
  typename LinearOperator::DomainFunctionType domainFunction( "domain Function", dSpace );
  typename LinearOperator::RangeFunctionType rangeFunction1( "range Function1", rSpace );
  typename LinearOperator::RangeFunctionType rangeFunction2( "range Function2", rSpace );

  // first test initialization
  RangeStencil< DomainSpaceType, RangeSpaceType > stencil( dSpace, rSpace, range );
  linOp.reserve( stencil );
  linOp.clear();

  Dune::Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > temp( dSpace, rSpace );

  // test assamble
  for( const auto &entry : range )
  {
    auto domainEntity = entry.first;
    auto rangeEntity = entry.second;

    {
      // check {add,set,addScaled}LocalMatrix
      temp.init( domainEntity, rangeEntity );
      temp.clear();
      for( const auto &p : permutation )
        temp.set( p.first, p.second, 1.0 );

      linOp.setLocalMatrix( domainEntity, rangeEntity, temp );

      linOp.addLocalMatrix( domainEntity, rangeEntity, temp );
      linOp.addScaledLocalMatrix( domainEntity, rangeEntity, temp, -1.0 );
    }
  }

  // finalize assamble
  linOp.communicate();


  // test getLocalMatrix
  for( const auto &entry : range )
  {
    auto domainEntity = entry.first;
    auto rangeEntity = entry.second;
    temp.init( domainEntity, rangeEntity );
    temp.clear();
    linOp.getLocalMatrix( domainEntity, rangeEntity, temp );

    for( const auto &p : permutation )
      if( temp.get( p.first, p.second ) != 1.0 )
        DUNE_THROW( Dune::NotImplemented, "LocalMatrix not set correctly" );

    // check old localMatrix implementation
    LocalMatrixType localMatrix = linOp.localMatrix( domainEntity, rangeEntity );

    for( const auto &p : permutation )
    {
      double val = localMatrix.get( p.first, p.second );
      if( val != 1.0 )
        DUNE_THROW( Dune::NotImplemented, "LocalMatrix not set correctly" );
    }
  }


  domainFunction.clear();
  int minNumDofs = std::min( domainFunction.blocks(), rangeFunction1.blocks() );
  for( int i = 0; i < minNumDofs; ++i )
  {
    auto block = domainFunction.dofVector()[ i ];
    for( unsigned int j = 0; j < block.size(); ++j )
      block[ j ] = i;
  }

  // test application
  linOp( domainFunction, rangeFunction1 );

  // mimic operator apply manually
  rangeFunction2.clear();
  Dune::DynamicVector< double > tmp;
  std::vector< double > values;
  for( const auto &entry : range )
  {
    auto domainEntity = entry.first;
    auto rangeEntity = entry.second;

    std::size_t oldSize = domainFunction.localFunction( domainEntity ).size();
    tmp.resize( oldSize );
    domainFunction.getLocalDofs( domainEntity, tmp );

    values.clear();
    std::size_t newSize = rangeFunction2.localFunction( rangeEntity ).size();
    values.resize( newSize, 0.0 );

    for( const auto &p : permutation )
      values[ p.first ] = tmp[ p.second ];
    rangeFunction2.setLocalDofs( rangeEntity, values );
  }
  rangeFunction2.communicate();

  rangeFunction1.dofVector() -= rangeFunction2.dofVector();

  double diff = rangeFunction1.scalarProductDofs( rangeFunction1 );
  if( std::abs( diff ) > 1e-8 )
    DUNE_THROW( Dune::NotImplemented, "Id Operator not assembled correctly" );
}

// Main Program
// ------------

int main ( int argc, char **argv )
try
{
  // initialize MPI manager and PETSc
  Dune::Fem::MPIManager::initialize( argc, argv );

  GridType grid( {1, 1}, {{2, 2}} );

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

#if not USE_PETSC && HAVE_PETSC
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
#endif // #if not USE_PETSC

  return 0;
}
catch( const Dune::Exception &exception )
{
  // display the exception message on the console
  std::cerr << exception << std::endl;
  return 1;
}
