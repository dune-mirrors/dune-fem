#ifndef DUNE_FEM_OPERATOR_LINEAR_BLOCKDIAGONAL_HH
#define DUNE_FEM_OPERATOR_LINEAR_BLOCKDIAGONAL_HH

#include <cassert>

#include <dune/common/fmatrix.hh>
#include <dune/common/nullptr.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/matrix/columnobject.hh>
#include <dune/fem/storage/objectstack.hh>

namespace Dune
{

  namespace Fem
  {

    // BlockDiagonalLinearOperator
    // ---------------------------

    template< class DiscreteFunctionSpace >
    class BlockDiagonalLinearOperator
    : public Fem::AssembledOperator< AdaptiveDiscreteFunction< DiscreteFunctionSpace > >
    {
      typedef BlockDiagonalLinearOperator< DiscreteFunctionSpace > ThisType;
      typedef Fem::AssembledOperator< AdaptiveDiscreteFunction< DiscreteFunctionSpace > > BaseType;

    public:
      typedef typename BaseType::DomainFunctionType DomainFunctionType;
      typedef typename BaseType::RangeFunctionType RangeFunctionType;

      typedef typename BaseType::DomainFieldType DomainFieldType;
      typedef typename BaseType::RangeFieldType RangeFieldType;

      typedef typename DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename RangeFunctionType::DiscreteFunctionSpaceType RangeSpaceType;

      typedef typename DomainSpaceType::EntityType DomainEntityType;
      typedef typename RangeSpaceType::EntityType RangeEntityType;

    private:
      static const int localBlockSize = DomainSpaceType::localBlockSize;

      class LocalMatrixTraits;
      class LocalMatrix;
      struct LocalMatrixFactory;

    public:
      typedef ObjectStack< LocalMatrixFactory > LocalMatrixStackType;
      typedef LocalMatrixWrapper< LocalMatrixStackType > LocalMatrixType;

      typedef ColumnObject< ThisType > LocalColumnObjectType;

      typedef Dune::FieldMatrix< RangeFieldType, localBlockSize, localBlockSize > LocalBlockType;

      BlockDiagonalLinearOperator ( const std::string &name,
                                    const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace )
      : name_( name ),
        space_( domainSpace ),
        localMatrixFactory_( *this ),
        localMatrixStack_( localMatrixFactory_ )
      {
        if( &domainSpace != &rangeSpace )
          DUNE_THROW( InvalidStateException, "BlockDiagonalLinearOperator must be created with identical spaces." );
      }

      void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
      {
        typedef typename std::vector< LocalBlockType >::const_iterator Iterator;
        const RangeFieldType *uit = u.leakPointer();
        RangeFieldType *wit = w.leakPointer();
        for( Iterator dit = diagonal_.begin(); dit != diagonal_.end(); ++dit )
        {
          dit->mv( uit, wit );
          uit += dit->M();
          wit += dit->N();
        }
        assert( uit == u.leakPointer() + u.size() );
        assert( wit == w.leakPointer() + w.size() );
      }

      void clear ()
      {
        typedef typename std::vector< LocalBlockType >::iterator Iterator;
        for( Iterator dit = diagonal_.begin(); dit != diagonal_.end(); ++dit )
          *dit = RangeFieldType( 0 );
      }

      template< class Functor >
      void forEach ( const Functor &functor )
      {
        typedef typename std::vector< LocalBlockType >::iterator Iterator;
        for( Iterator dit = diagonal_.begin(); dit != diagonal_.end(); ++dit )
          functor( *dit );
      }

      void invert ()
      {
        typedef typename std::vector< LocalBlockType >::iterator Iterator;
        for( Iterator dit = diagonal_.begin(); dit != diagonal_.end(); ++dit )
          dit->invert();
      }

      void communicate () {}

      template< class Stencil >
      void reserve ( const Stencil &stencil, bool verbose = false )
      {
        // note: to use DynamicMatrix as LocalBlockType, also resize local blocks
        diagonal_.resize( domainSpace().blockMapper().size() );
      }

      LocalColumnObjectType localColumn ( const DomainEntityType &domainEntity ) const
      {
        return LocalColumnObjectType( *this, domainEntity );
      }

      LocalMatrixType localMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const;

      const DomainSpaceType &domainSpace () const { return space_; }
      const RangeSpaceType &rangeSpace () const { return space_; }

      const std::string &name () const { return name_; }

    private:
      std::string name_;
      const RangeSpaceType &space_;
      std::vector< LocalBlockType > diagonal_;

      LocalMatrixFactory localMatrixFactory_;
      mutable LocalMatrixStackType localMatrixStack_;
    };



    // BlockDiagonalLinearOperator::LocalMatrixTraits
    // ----------------------------------------------

    template< class DiscreteFunctionSpace >
    class BlockDiagonalLinearOperator< DiscreteFunctionSpace >::LocalMatrixTraits
    {
      typedef BlockDiagonalLinearOperator< DiscreteFunctionSpace > OperatorType;

    public:
      typedef typename OperatorType::LocalMatrix LocalMatrixType;

      typedef typename OperatorType::RangeFieldType RangeFieldType;

      typedef typename OperatorType::DomainSpaceType DomainSpaceType;
      typedef typename OperatorType::RangeSpaceType RangeSpaceType;

      typedef RangeFieldType LittleBlockType;
    };



    // BlockDiagonalLinearOperator::LocalMatrix
    // ----------------------------------------

    template< class DiscreteFunctionSpace >
    class BlockDiagonalLinearOperator< DiscreteFunctionSpace >::LocalMatrix
    : public LocalMatrixInterface< LocalMatrixTraits >
    {
      typedef LocalMatrix ThisType;
      typedef LocalMatrixInterface< LocalMatrixTraits > BaseType;

    public:
      typedef BlockDiagonalLinearOperator< DiscreteFunctionSpace > OperatorType;

      typedef typename BaseType::RangeFieldType RangeFieldType;

      typedef typename BaseType::DomainBasisFunctionSetType DomainBasisFunctionSetType;
      typedef typename BaseType::RangeBasisFunctionSetType RangeBasisFunctionSetType;

      typedef typename BaseType::DomainEntityType DomainEntityType;
      typedef typename BaseType::RangeEntityType RangeEntityType;

    private:
      typedef DomainBasisFunctionSetType BasisFunctionSetType;
      typedef typename OperatorType::LocalBlockType LocalBlockType;

      struct SetLocalBlockFunctor
      {
        SetLocalBlockFunctor ( std::vector< LocalBlockType > &diagonal,
                               LocalBlockType *&localBlock )
        : diagonal_( diagonal ),
          localBlock_( localBlock )
        {}

        void operator() ( int localDoF, std::size_t globalDoF )
        {
          assert( localDoF == 0 );
          assert( globalDoF < diagonal_.size() );
          localBlock_ = &diagonal_[ globalDoF ];
        }

      private:
        std::vector< LocalBlockType > &diagonal_;
        LocalBlockType *&localBlock_;
      };

    public:
      explicit LocalMatrix ( OperatorType &op )
      : op_( &op ),
        localBlock_( nullptr )
      {}

      void init ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        basisFunctionSet_ = domainSpace().basisFunctionSet( domainEntity );
        SetLocalBlockFunctor f( op_->diagonal_, localBlock_ );
        domainSpace().blockMapper().mapEach( domainEntity, f );
      }

      void clear () { localBlock() = RangeFieldType( 0 ); }
      void scale ( const RangeFieldType &a ) { localBlock() *= a; }

      RangeFieldType get ( int i, int j ) const { return localBlock()[ i ][ j ]; }
      void add ( int i, int j, const RangeFieldType &value ) { localBlock()[ i ][ j ] += value; }
      void set ( int i, int j, const RangeFieldType &value ) { localBlock()[ i ][ j ] = value; }

      void clearRow ( int i ) { localBlock()[ i ] = RangeFieldType( 0 ); }

      void clearCol ( int j )
      {
        for( int i = 0; i < rows(); ++i )
          localBlock()[ i ][ j ] = RangeFieldType( 0 );
      }

      void unitRow ( int i ) DUNE_DEPRECATED
      {
        localBlock()[ i ] = RangeFieldType( 0 );
        localBlock()[ i ][ i ] = RangeFieldType( 1 );
      }

      template< class DomainLocalFunction, class RangeLocalFunction >
      void multiplyAdd ( const DomainLocalFunction &x, RangeLocalFunction &y ) const
      {
        localBlock().umv( x, y );
      }

      void finalize () {}
      void resort () {}

      int rows () const { return localBlock().N(); }
      int columns () const { return localBlock().M(); }

      const DomainSpaceType &domainSpace () const { return op_->domainSpace(); }
      const RangeSpaceType &rangeSpace () const { return op_->rangeSpace(); }

      const DomainBasisFunctionSetType &domainBasisFunctionSet () const { return basisFunctionSet_; }
      const RangeBasisFunctionSetType &rangeBasisFunctionSet () const { return basisFunctionSet_; }

      const DomainEntityType &domainEntity () const { return domainBasisFunctionSet().entity(); }
      const RangeEntityType &rangeEntity () const { return rangeBasisFunctionSet().entity(); }

    private:
      const LocalBlockType &localBlock () const { assert( localBlock_ ); return *localBlock_; }
      LocalBlockType &localBlock () { assert( localBlock_ ); return *localBlock_; }

      OperatorType *op_;
      BasisFunctionSetType basisFunctionSet_;
      LocalBlockType *localBlock_;
    };



    // BlockDiagonalLinearOperator::LocalMatrixFactory
    // -----------------------------------------------

    template< class DiscreteFunctionSpace >
    struct BlockDiagonalLinearOperator< DiscreteFunctionSpace >::LocalMatrixFactory
    {
      typedef BlockDiagonalLinearOperator< DiscreteFunctionSpace > OperatorType;
      typedef LocalMatrix ObjectType;

      explicit LocalMatrixFactory ( OperatorType &op )
      : op_( &op )
      {}

      ObjectType *newObject () const { return new ObjectType( *op_ ); }

    private:
      OperatorType *op_;
    };



    // Implementation of BlockDiagonalLinearOperator
    // ---------------------------------------------

    template< class DiscreteFunctionSpace >
    inline typename BlockDiagonalLinearOperator< DiscreteFunctionSpace >::LocalMatrixType
    BlockDiagonalLinearOperator< DiscreteFunctionSpace >
      ::localMatrix ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity ) const
    {
      return LocalMatrixType( localMatrixStack_, domainEntity, rangeEntity );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_LINEAR_BLOCKDIAGONAL_HH
