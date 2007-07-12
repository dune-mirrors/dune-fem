#ifndef DUNE_FEM_INTEGRATIONOPERATOR_HH
#define DUNE_FEM_INTEGRATIONOPERATOR_HH

#include <dune/fem/function/common/temporarylocalfunction.hh>
#include <dune/fem/operator/matrix/localmatrix.hh>

namespace Dune
{

  template< class LocalOperatorImp, class RangeDiscreteFunctionImp >
  class DefaultIntegrationOperatorTraits
  {
  public:
    typedef LocalOperatorImp LocalOperatorType;

    typedef RangeDiscreteFunctionImp RangeFunctionType;

  private:
    typedef DefaultIntegrationOperatorTraits
      < LocalOperatorType, RangeFunctionType >
      ThisType;

  public:
    typedef typename LocalOperatorType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename LocalOperatorType :: DomainFunctionType DomainFunctionType;
    
    typedef typename RangeFunctionType :: FunctionSpaceType RangeFunctionSpaceType;
    
    typedef typename DomainFunctionSpaceType :: RangeFieldType DomainFieldType;
    typedef typename RangeFunctionSpaceType :: RangeFieldType RangeFieldType;
  };



  template< class TraitsImp, bool isLinear >
  class IntegrationOperator;



  template< class TraitsImp >
  class IntegrationOperator< TraitsImp, true >
  : public Operator< typename TraitsImp :: DomainFieldType,
                     typename TraitsImp :: RangeFieldType,
                     typename TraitsImp :: DomainFunctionType,
                     typename TraitsImp :: RangeFunctionType >
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef IntegrationOperator< TraitsType, true > ThisType;

  public:
    typedef typename TraitsType :: LocalOperatorType LocalOperatorType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename LocalOperatorType :: DomainFunctionType
      DomainFunctionType;
    typedef typename LocalOperatorType :: DomainFunctionSpaceType
      DomainFunctionSpaceType;

    typedef typename RangeFunctionType :: FunctionSpaceType
      RangeFunctionSpaceType;

    typedef TemporaryLocalFunction< RangeFunctionSpaceType >
      RangeTemporaryLocalFunctionType;
    typedef TemporaryLocalMatrix< DomainFunctionSpaceType, RangeFunctionSpaceType >
      TemporaryLocalMatrixType;

  protected:
    const LocalOperatorType &localOperator_;
    const DomainFunctionSpaceType &domainFunctionSpace_;
    const RangeFunctionSpaceType &rangeFunctionSpace_;
    
  public:
    inline IntegrationOperator ( const LocalOperatorType &localOperator,
                                 const DomainFunctionSpaceType &domainFunctionSpace,
                                 const RangeFunctionSpaceType &rangeFunctionSpace )
    : localOperator_( localOperator ),
      domainFunctionSpace_( domainFunctionSpace ),
      rangeFunctionSpace_( rangeFunctionSpace )
    {
    }

    inline void operator()
      ( const DomainFunctionType &u,
        RangeFunctionType &w
      ) const
    {
      // type of iterator over grid partition for range function space
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;
     
      // type of base function set for range function space
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

      // type of local function for range functions
      typedef typename RangeFunctionType :: LocalFunctionType
        RangeLocalFunctionType;

      // obtain range function space (discrete function space)
      const RangeFunctionSpaceType &rangeFunctionSpace = w.space();

      // create a temporary local function
      RangeTemporaryLocalFunctionType w_temp( rangeFunctionSpace );
      
      // clear destination function
      w.clear();
     
      // Loop over entire grid
      const IteratorType &end = rangeFunctionSpace.end();
      for( IteratorType it = rangeFunctionSpace.begin(); it != end; ++it )
      {
        // apply the local operator
        w_temp.init( *it );
        localOperator_( *it, u, w_temp );

        // update the destination function
        RangeLocalFunctionType w_local = w.localFunction( *it );
        const int numDofs = w_local.numDofs();
        for( int i = 0; i < numDofs; ++i )
          w_local[ i ] += w_temp[ i ];
      }
    }

    template< class MatrixType >
    inline void assembleMatrix ( MatrixType &matrix ) const
    {
      // type of iterator over grid partition
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;
     
      // type of base function sets
      typedef typename DomainFunctionSpaceType :: BaseFunctionSetType
        DomainBaseFunctionSetType;
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        RangeBaseFunctionSetType;

      // row type for local matrix
      typedef typename TemporaryLocalMatrixType :: RowType LocalRowType;

      const DomainFunctionSpaceType &domainFunctionSpace = this->domainFunctionSpace();
      const RangeFunctionSpaceType &rangeFunctionSpace = this->rangeFunctionSpace();

      // create local matrix
      TemporaryLocalMatrixType localMatrix( domainFunctionSpace, rangeFunctionSpace );
      
      // clear global matrix
      matrix.clear();

      // Loop over entire grid
      const IteratorType &end = rangeFunctionSpace.end();
      for( IteratorType it = rangeFunctionSpace.begin(); it != end; ++it )
      {
        // obtain local matrix from local operator
        localMatrix.init( *it );
        localOperator_.assembleMatrix( *it, localMatrix );

        // update global matrix
        const unsigned int numLocalRows = localMatrix.rows();
        const unsigned int numLocalColumns = localMatrix.columns();
        for( unsigned int i = 0; i < numLocalRows; ++i ) {
          const LocalRowType localRow = localMatrix[ i ];

          const int globalRowIndex = rangeFunctionSpace.mapToGlobal( *it, i );
          for( unsigned int j = 0; j < numLocalColumns; ++j ) {
            const int globalColIndex = domainFunctionSpace.mapToGlobal( *it, j );
            matrix.add( globalRowIndex, globalColIndex, localRow[ j ] );
          }
        }
      }
    }

    inline const DomainFunctionSpaceType &domainFunctionSpace () const
    {
      return domainFunctionSpace_;
    }

    inline const RangeFunctionSpaceType &rangeFunctionSpace () const
    {
      return rangeFunctionSpace_;
    }
  };


  
  template< class TraitsImp >
  class IntegrationOperator< TraitsImp, false >
  : public Operator< typename TraitsImp :: DomainFieldType,
                     typename TraitsImp :: RangeFieldType,
                     typename TraitsImp :: DomainFunctionType,
                     typename TraitsImp :: RangeFunctionType >
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef IntegrationOperator< TraitsType, false > ThisType;

  public:
    typedef typename TraitsType :: LocalOperatorType LocalOperatorType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename LocalOperatorType :: DomainFunctionType
      DomainFunctionType;
    typedef typename LocalOperatorType :: DomainFunctionSpaceType
      DomainFunctionSpaceType;

    typedef typename RangeFunctionType :: FunctionSpaceType
      RangeFunctionSpaceType;

    typedef TemporaryLocalFunction< RangeFunctionSpaceType >
      RangeTemporaryLocalFunctionType;
    typedef TemporaryLocalMatrix< DomainFunctionSpaceType, RangeFunctionSpaceType >
      TemporaryLocalMatrixType;

  protected:
    const LocalOperatorType &localOperator_;
    
  public:
    inline IntegrationOperator ( const LocalOperatorType &localOperator )
    : localOperator_( localOperator )
    {
    }

    inline void operator()
      ( const DomainFunctionType &u,
        RangeFunctionType &w
      ) const
    {
      // type of iterator over grid partition for range function space
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;
     
      // type of base function set for range function space
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

      // type of local function for range functions
      typedef typename RangeFunctionType :: LocalFunctionType
        RangeLocalFunctionType;

      // obtain range function space (discrete function space)
      const RangeFunctionSpaceType &rangeFunctionSpace = w.space();

      // create a temporary local function
      RangeTemporaryLocalFunctionType w_temp( rangeFunctionSpace );
      
      // clear destination function
      w.clear();
     
      // Loop over entire grid
      const IteratorType &end = rangeFunctionSpace.end();
      for( IteratorType it = rangeFunctionSpace.begin(); it != end; ++it )
      {
        // apply the local operator
        w_temp.init( *it );
        localOperator_( *it, u, w_temp );

        // update the destination function
        RangeLocalFunctionType w_local = w.localFunction( *it );
        const int numDofs = w_local.numDofs();
        for( int i = 0; i < numDofs; ++i )
          w_local[ i ] += w_temp[ i ];
      }
    }
  };

  
}

#endif
