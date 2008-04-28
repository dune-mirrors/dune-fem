#ifndef DUNE_CACHEDLINEAROPERATOR_HH
#define DUNE_CACHEDLINEAROPERATOR_HH

#include <dune/fem/operator/common/operator.hh>

namespace Dune
{
  //! @ingroup LinearDiscOperatos 
  //! \newimplementation see fem/operator/concept.tex
  template< class WrappedOperatorImp, class SystemMatrixImp >
  class CachedLinearOperator
  : public Operator< typename WrappedOperatorImp :: DomainFieldType,
                     typename WrappedOperatorImp :: RangeFieldType,
                     typename WrappedOperatorImp :: DomainFunctionType,
                     typename WrappedOperatorImp :: RangeFunctionType >
  {
  public:
    typedef WrappedOperatorImp WrappedOperatorType;
    typedef SystemMatrixImp SystemMatrixType;

  private:
    typedef CachedLinearOperator< WrappedOperatorType, SystemMatrixType > ThisType;

  public:
    typedef typename WrappedOperatorType :: DomainFieldType DomainFieldType;
    typedef typename WrappedOperatorType :: RangeFieldType RangeFieldType;
   
    typedef typename WrappedOperatorType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename WrappedOperatorType :: RangeFunctionSpaceType RangeFunctionSpaceType;
    
    typedef typename WrappedOperatorType :: DomainFunctionType DomainFunctionType;
    typedef typename WrappedOperatorType :: RangeFunctionType RangeFunctionType;

    typedef typename WrappedOperatorType :: DomainProjectionType DomainProjectionType;
    typedef typename WrappedOperatorType :: RangeProjectionType RangeProjectionType;

  protected:
    const WrappedOperatorType *wrappedOperator_;
    mutable SystemMatrixType *matrix_;

  public:
    inline CachedLinearOperator ()
    : wrappedOperator_( NULL ),
      matrix_( NULL )
    {
    }

    inline CachedLinearOperator ( const WrappedOperatorType &wrappedOperator )
    : wrappedOperator_( &wrappedOperator ),
      matrix_( NULL )
    {
    }

    inline ~CachedLinearOperator ()
    {
      if( matrix_ != NULL )
        delete matrix_;
    }

    inline ThisType &operator= ( const WrappedOperatorType &wrappedOperator )
    {
      wrappedOperator_ = &wrappedOperator;
      rebuild();
    }

    inline void operator() ( const DomainFunctionType &u, RangeFunctionType &w ) const
    {
      assemble();
      matrix_->apply( u, w ); 
    }

    template< class MatrixType >
    inline void assembleMatrix ( MatrixType &matrix ) const
    {
      assemble();
      DUNE_THROW( NotImplemented,
                  "CachedLinearOperator :: assembleMatrix is not implemented, yet." );
    }

    inline void rebuild ()
    {
      if( matrix_ != NULL )
      {
        delete matrix_;
        matrix_ = NULL;
      }
    }

    inline const DomainFunctionSpaceType &domainFunctionSpace () const
    {
      return wrappedOperator_->domainFunctionSpace();
    }

    inline const RangeFunctionSpaceType &rangeFunctionSpace () const
    {
      return wrappedOperator_->rangeFunctionSpace();
    }
    
    inline const DomainProjectionType domainProjection () const
    {
      return wrappedOperator_->domainProjection();
    }

    inline const RangeProjectionType rangeProjection () const
    {
      return wrappedOperator_->rangeProjection();
    }
    
    inline const SystemMatrixType& systemMatrix () const
    {
      assemble();
      return *matrix_;
    }

  private:
    inline void assemble () const
    {
      if( matrix_ != NULL )
        return;

      const DomainFunctionSpaceType &domainFunctionSpace
        = this->domainFunctionSpace();
      const RangeFunctionSpaceType &rangeFunctionSpace
        = this->rangeFunctionSpace();

      const unsigned int rows = rangeFunctionSpace.size();
      const unsigned int columns = domainFunctionSpace.size();
      matrix_ = new SystemMatrixType( rows, columns, 50 );
      assert( matrix_ != NULL );

      wrappedOperator_->assembleMatrix( *matrix_ );
    }
  };
  
}

#endif
