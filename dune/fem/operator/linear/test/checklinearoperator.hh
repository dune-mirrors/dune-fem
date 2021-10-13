#ifndef DUNE_FEM_OPERATOR_LINEAR_TEST_CHECKLINEAROPERATOR_HH
#define DUNE_FEM_OPERATOR_LINEAR_TEST_CHECKLINEAROPERATOR_HH

#include <cstddef>
#include <cmath>

#include <algorithm>
#include <utility>
#include <vector>

#include <dune/common/classname.hh>
#include <dune/common/dynvector.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/typeutilities.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/operator/common/stencil.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>

namespace Dune
{

  namespace Fem
  {

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

        Iterator () = default;

        Iterator ( Iterator1 it1, Iterator2 it2 ) : it_( it1, it2 ) {}

        Iterator &operator++ () { ++it_.first; ++it_.second; return *this; }

        std::pair< Entity1, Entity2 > operator* () const { return std::make_pair( *it_.first, *it_.second ); }

        bool operator!= ( const Iterator &other ) const { return (it_ != other.it_); }
        bool operator== ( const Iterator &other ) const { return (it_ == other.it_); }

      private:
        std::pair< Iterator1, Iterator2 > it_;
      };

      DiagonalRange ( const DomainSpace &dSpace, const RangeSpace &rSpace )
        : dSpace_( dSpace ), rSpace_( rSpace )
      {}

      Iterator begin () const { return Iterator( dSpace_.begin(), rSpace_.begin() ); }
      Iterator end () const { return Iterator( dSpace_.end(), rSpace_.end() ); }

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

      virtual void setupStencil() const
      {
        DUNE_THROW(NotImplemented,"This is done in the constructor");
      }
    };



    namespace CheckLinearOperator
    {

      template< class LocalMatrix >
      inline void verifyPermutation ( const LocalMatrix &localMatrix, const std::vector< std::pair< int, int > > &permutation )
      {
        for( const auto &p : permutation )
        {
          if( localMatrix.get( p.first, p.second ) != 1.0 )
            DUNE_THROW( Dune::NotImplemented, "LocalMatrix not set correctly" );
        }
      }

      template< class LinearOperator, class DomainEntity, class RangeEntity >
      inline void_t< typename LinearOperator::LocalMatrixType >
      verifyLocalMatrixPermutation ( const LinearOperator &linOp, const DomainEntity &domainEntity, const RangeEntity &rangeEntity,
                                     const std::vector< std::pair< int, int > > &permutation, PriorityTag< 1 > )
      {
        Dune::Fem::TemporaryLocalMatrix< typename LinearOperator::DomainSpaceType, typename LinearOperator::RangeSpaceType >
          localMat( linOp.domainSpace(), linOp.rangeSpace() );
        localMat.bind( domainEntity, rangeEntity );
        linOp.getLocalMatrix( domainEntity, rangeEntity, localMat );
        verifyPermutation( localMat, permutation );
        localMat.unbind();
      }

      template< class LinearOperator, class DomainEntity, class RangeEntity >
      inline void verifyLocalMatrixPermutation ( const LinearOperator &linOp, const DomainEntity &domainEntity, const RangeEntity &rangeEntity,
                                                 const std::vector< std::pair< int, int > > &permutation, PriorityTag< 0 > )
      {
        static bool warn = true;
        if( warn )
          std::cout << "Note: " << className( linOp ) << " does not support old-style local matrix interface." << std::endl;
        warn = false;
      }

      template< class LinearOperator, class DomainEntity, class RangeEntity >
      inline void verifyLocalMatrixPermutation ( const LinearOperator &linOp, const DomainEntity &domainEntity, const RangeEntity &rangeEntity,
                                                const std::vector< std::pair< int, int > > &permutation )
      {
        verifyLocalMatrixPermutation( linOp, domainEntity, rangeEntity, permutation, PriorityTag< 42 >() );
      }

    } // namespace CheckLinearOperator



    // checkLinearOperator
    // -------------------

    template< class LinearOperator, class Range >
    inline void checkLinearOperator ( LinearOperator &linOp, const Range &range, const std::vector< std::pair< int, int > > &permutation )
    {
      typedef typename LinearOperator::DomainFunctionType DomainFunctionType;
      typedef typename LinearOperator::RangeFunctionType RangeFunctionType;

      // get type of domain and range function
      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      // get spaces
      const DomainSpaceType &dSpace = linOp.domainSpace();
      const RangeSpaceType &rSpace = linOp.rangeSpace();

      // construct some domain and range functions
      DomainFunctionType domainFunction( "domain Function", dSpace );
      RangeFunctionType rangeFunction1( "range Function1", rSpace );
      RangeFunctionType rangeFunction2( "range Function2", rSpace );

      // first test initialization
      RangeStencil< DomainSpaceType, RangeSpaceType > stencil( dSpace, rSpace, range );

      for( int k=0; k<4; ++k )
      {
        // check setting operator entries to zero
        linOp.reserve( stencil );
        linOp.clear();

        Dune::Fem::TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType > temp( dSpace, rSpace );

        // test assemble
        // check {add,addScaled}LocalMatrix
        for( const auto &entry : range )
        {
          auto domainEntity = entry.first;
          auto rangeEntity = entry.second;

          {
            temp.init( domainEntity, rangeEntity );
            temp.clear();
            for( const auto &p : permutation )
              temp.set( p.first, p.second, 1.0 );

            linOp.addLocalMatrix( domainEntity, rangeEntity, temp );
            linOp.addScaledLocalMatrix( domainEntity, rangeEntity, temp, -1.0 );
          }
        }

        linOp.flushAssembly();

        // check {set}LocalMatrix
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
          }
        }

        linOp.flushAssembly();
        if (k>=2) // just for checking 'clearRow' and reassembly so don't check content of matrix
        {
          try
          {
            linOp.unitRow(0);
            linOp.finalize();
          }
          catch( const Dune::NotImplemented &exception )
          {
            std::cout << "operator does not implement unitRow method\n";
          }
          continue;
        }

        // test getLocalMatrix
        for( const auto &entry : range )
        {
          auto domainEntity = entry.first;
          auto rangeEntity = entry.second;

          temp.init( domainEntity, rangeEntity );
          temp.clear();
          linOp.getLocalMatrix( domainEntity, rangeEntity, temp );
          CheckLinearOperator::verifyPermutation( temp, permutation );

          // check old localMatrix implementation
          CheckLinearOperator::verifyLocalMatrixPermutation( linOp, domainEntity, rangeEntity, permutation );
        }

        // test application
        domainFunction.clear();
        std::size_t minNumDofs = std::min( domainFunction.blocks(), rangeFunction1.blocks() );
        for( std::size_t i = 0; i < minNumDofs; ++i )
        {
          auto block = domainFunction.dofVector()[ i ];
          Hybrid::forEach( typename DomainSpaceType::LocalBlockIndices(), [ &block, i ] ( auto &&j ) { block[ j ] = i; } );
        }
        linOp( domainFunction, rangeFunction1 );

        // mimic operator apply manually
        rangeFunction2.clear();
        Dune::DynamicVector< double > tmp;
        std::vector< double > values;
        for( const auto &entry : range )
        {
          auto domainEntity = entry.first;
          auto rangeEntity = entry.second;

          tmp.resize( dSpace.blockMapper().numDofs( domainEntity ) * DomainSpaceType::localBlockSize );
          domainFunction.getLocalDofs( domainEntity, tmp );

          values.clear();
          values.resize( rSpace.blockMapper().numDofs( rangeEntity ) * RangeSpaceType::localBlockSize, 0.0 );

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
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_TEST_CHECKLINEAROPERATOR_HH
