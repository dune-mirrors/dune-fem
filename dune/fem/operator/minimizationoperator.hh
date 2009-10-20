#ifndef DUNE_FEM_MINIMIZATIONOPERATOR_HH
#define DUNE_FEM_MINIMIZATIONOPERATOR_HH

#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/function/localfunction/temporarylocalfunction.hh>

namespace Dune
{

  template< class LocalOperatorImp, class RangeDiscreteFunctionImp >
  class DefaultMinimizationOperatorTraits
  {
  public:
    typedef LocalOperatorImp LocalOperatorType;

    typedef RangeDiscreteFunctionImp RangeFunctionType;

  private:
    typedef DefaultMinimizationOperatorTraits
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



  template< class TraitsImp >
  class MinimizationOperator
  : public Operator< typename TraitsImp :: DomainFieldType,
                     typename TraitsImp :: RangeFieldType,
                     typename TraitsImp :: DomainFunctionType,
                     typename TraitsImp :: RangeFunctionType >
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef MinimizationOperator< TraitsType > ThisType;

  public:
    typedef typename TraitsType :: LocalOperatorType LocalOperatorType;

    typedef typename TraitsType :: DomainFunctionType DomainFunctionType;
    typedef typename TraitsType :: RangeFunctionType RangeFunctionType;

    typedef typename TraitsType :: DomainFunctionSpaceType DomainFunctionSpaceType;
    typedef typename TraitsType :: RangeFunctionSpaceType RangeFunctionSpaceType;

    typedef typename TraitsType :: DomainFieldType DomainFieldType;
    typedef typename TraitsType :: RangeFieldType RangeFieldType;

  protected:
    typedef typename RangeFunctionType :: LocalFunctionType RangeLocalFunctionType;
    
  protected:
    const LocalOperatorType &localOperator_;

  public:
    inline explicit MinimizationOperator ( const LocalOperatorType &localOperator )
    : localOperator_( localOperator )
    {
    }
    
    inline void operator () ( const DomainFunctionType &u,
                              RangeFunctionType &w ) const
    {
      typedef typename RangeFunctionSpaceType :: IteratorType IteratorType;
      typedef typename IteratorType :: Entity EntityType;

      typedef typename RangeFunctionSpaceType :: LocalFunctionType
        TemporaryLocalFunctionType;

      typedef typename RangeFunctionType :: LocalFunctionType RangeLocalFunctionType;
      typedef typename RangeFunctionType :: DofIteratorType DofIteratorType;

      const RangeFunctionSpaceType &rangeFunctionSpace = w.space();

      const DofIteratorType dend = w.dend();
      for( DofIteratorType dit = w.dbegin(); dit != dend; ++dit )
        (*dit) = (RangeFieldType)INFINITY;

      const IteratorType end = rangeFunctionSpace.end();
      for( IteratorType it = rangeFunctionSpace.begin(); it != end; ++it )
      {
        const EntityType &entity = *it;

        TemporaryLocalFunctionType w_temp = rangeFunctionSpace.localFunction( entity );
        localOperator_( entity, u, w_temp );

        RangeLocalFunctionType w_local = w.localFunction( entity );
        const unsigned int numDofs = w_local.numDofs();
        for( unsigned int i = 0; i < numDofs; ++i )
        {
          RangeFieldType &dof = w_local[ i ];
          RangeFieldType &update = w_temp[ i ];
          dof = dof < update ? dof : update;
        }
      }
    }
  };
  
}

#endif
