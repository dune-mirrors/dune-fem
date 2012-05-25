#ifndef DUNE_FEM_INTEGRATIONOPERATOR_HH
#define DUNE_FEM_INTEGRATIONOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>

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
    typedef typename LocalOperatorType :: DomainFunctionSpaceType
      DomainFunctionSpaceType;
    typedef typename LocalOperatorType :: DomainFunctionType DomainFunctionType;
    
    typedef typename RangeFunctionType :: DiscreteFunctionSpaceType
      RangeFunctionSpaceType;
    
    typedef typename LocalOperatorType :: DomainFieldType DomainFieldType;
    typedef typename LocalOperatorType :: RangeFieldType RangeFieldType;
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

    typedef typename TraitsType :: DomainFunctionType DomainFunctionType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename TraitsType :: DomainFunctionSpaceType
      DomainFunctionSpaceType;
    typedef typename TraitsType :: RangeFunctionSpaceType
      RangeFunctionSpaceType;

    typedef typename TraitsType :: DomainFieldType DomainFieldType;
    typedef typename TraitsType :: RangeFieldType RangeFieldType;

  protected:
    typedef typename RangeFunctionType :: LocalFunctionType
      RangeLocalFunctionType;

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
      typedef typename IteratorType :: Entity EntityType;
     
      // type of base function set for range function space
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

      // type of temporary local function (within range space)
      typedef typename RangeFunctionSpaceType :: LocalFunctionType
        TemporaryLocalFunctionType;

      // obtain range function space (discrete function space)
      const RangeFunctionSpaceType &rangeFunctionSpace = w.space();
     
      // clear destination function
      w.clear();
     
      // Loop over entire grid
      const IteratorType &end = rangeFunctionSpace.end();
      for( IteratorType it = rangeFunctionSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;
        
         // obtain a temporary local function
        TemporaryLocalFunctionType w_temp = rangeFunctionSpace.localFunction( entity );
          
        // apply the local operator
        localOperator_( entity, u, w_temp );

        // update the destination function
        RangeLocalFunctionType w_local = w.localFunction( entity );
        w_local += w_temp;
      }
    }

    template< class MatrixType >
    inline void assembleMatrix ( MatrixType &matrix ) const
    {
      // type of iterator over grid partition
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;

      typedef typename IteratorType :: Entity EntityType;
     
      // type of base function sets
      typedef typename DomainFunctionSpaceType :: BaseFunctionSetType
        DomainBaseFunctionSetType;
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        RangeBaseFunctionSetType;

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
        const EntityType &entity = *it;
        // obtain local matrix from local operator
        localMatrix.init( entity, entity );
        localOperator_.assembleMatrix( entity, localMatrix );

        // update global matrix
        const unsigned int numLocalRows = localMatrix.rows();
        const unsigned int numLocalColumns = localMatrix.columns();
        for( unsigned int i = 0; i < numLocalRows; ++i ) {
          const int globalRowIndex = rangeFunctionSpace.mapToGlobal( entity, i );
          for( unsigned int j = 0; j < numLocalColumns; ++j ) {
            const int globalColIndex = domainFunctionSpace.mapToGlobal( entity, j );
            matrix.add( globalRowIndex, globalColIndex, localMatrix.get( i, j ) );
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

    typedef typename TraitsType :: DomainFunctionType DomainFunctionType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename TraitsType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename TraitsType :: RangeFunctionSpaceType RangeFunctionSpaceType;

    typedef typename TraitsType :: DomainFieldType DomainFieldType;
    typedef typename TraitsType :: RangeFieldType RangeFieldType;

  protected:
    typedef typename RangeFunctionType :: LocalFunctionType
      RangeLocalFunctionType;
    
    typedef TemporaryLocalFunction< RangeFunctionSpaceType >
      RangeTemporaryLocalFunctionType;

  protected:
    const LocalOperatorType localOperator_;
    
  public:
    inline explicit IntegrationOperator ( const LocalOperatorType &localOperator )
    : localOperator_( localOperator )
    {
    }

    inline void operator () ( const DomainFunctionType &u,
                              RangeFunctionType &w ) const
    {
      // type of iterator over grid partition for range function space
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;
     
      // type of base function set for range function space
      typedef typename RangeFunctionSpaceType :: BaseFunctionSetType
        BaseFunctionSetType;

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
        w_local += w_temp;
      }
    }
  };
  
}

#endif
