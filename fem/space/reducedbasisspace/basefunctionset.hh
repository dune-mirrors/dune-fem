#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH

#include <dune/fem/space/common/basefunctioninterface.hh>

namespace Dune
{

  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSet;


  
  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSetTraits
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

  private:
    typedef ReducedBasisBaseFunctionSetTraits< BaseFunctionType > ThisType;

  public:
    typedef ReducedBasisBaseFunctionSet< BaseFunctionImp > BaseFunctionSetType;
    
    typedef typename BaseFunctionType :: FunctionSpaceType
      BaseFunctionSpaceType;
     
    typedef typename BaseFunctionSpaceType :: BaseFunctionSetType :: FunctionSpaceType
      FunctionSpaceType;
  };



  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSet
  : BaseFunctionSetDefault< ReducedBasisBaseFunctionSetTraits< BaseFunctionImp > >
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

    typedef ReducedBasisBaseFunctionSetTraits< BaseFunctionType > TraitsType;

  private:
    typedef ReducedBasisBaseFunctionSet< BaseFunctionType > ThisType;
    typedef BaseFunctionSetDefault< TraitsType > BaseType;

  public:
    typedef ThisType BaseFunctionSetType;
    
    typedef typename BaseFunctionType :: FunctionSpaceType
      BaseFunctionSpaceType;
    typedef typename BaseFunctionType :: LocalFunctionType
      LocalBaseFunctionType;

    typedef typename BaseFunctionSpaceType :: GridPartType GridPartType;
   
    typedef typename TraitsType :: FunctionSpaceType FunctionSpaceType;
    
    typedef typename FunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: RangeType RangeType;

    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType :: HessianRangeType HessianRangeType;

    enum { DimDomain = FunctionSpaceType :: DimDomain };
    enum { DimRange = FunctionSpaceType :: DimRange };

    typedef typename GridPartType :: EntityCodim0Type EntityCodim0Type;

    typedef std :: vector< BaseFunctionType* > BaseFunctionListType;

  protected:
    const BaseFunctionListType *baseFunctionList_;
    int numBaseFunctions_;
    
    const EntityCodim0Type *entity_;

  public:
    inline ReducedBasisBaseFunctionSet ()
    : baseFunctionList_( NULL ),
      entity_( NULL )
    {
    }

    inline ReducedBasisBaseFunctionSet ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( &baseFunctionList ),
      entity_( NULL )
    {
    }

    inline ReducedBasisBaseFunctionSet ( const BaseFunctionListType &baseFunctionList,
                                             const EntityCodim0Type &entity )
    : baseFunctionList_( &baseFunctionList ),
      entity_( &entity )
    {
    }

    inline ReducedBasisBaseFunctionSet ( const ThisType &other )
    : baseFunctionList_( other.baseFunctionList_ ),
      entity_( other.entity_ )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      baseFunctionList_ = other.baseFunctionList_;
      entity_ = other.entity_;
    }

    template< int diffOrd >
    inline void evaluate ( int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

      typedef typename LocalBaseFunctionType :: BaseFunctionSetType LocalBaseFunctionSetType;

      LocalBaseFunctionType localBaseFunction = baseFunctionList_[ baseFunction ]->localFunction( entity_ );
      LocalBaseFunctionSetType &localBaseFunctionSet = localBaseFunction.baseFunctionSet();
      const int numLocalBaseFunctions = localBaseFunctionSet.numBaseFunctions();
      
      phi = 0;
      for( int i = 0; i < numLocalBaseFunctions; ++i )
      {
        RangeType psi;
        localBaseFunctionSet.evaluate( i, diffVariable, x, psi );
        phi.axpy( localBaseFunction[ i ], psi );
      }
    }

    template< int diffOrd, class QuadratureType >
    inline void evaluate ( int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           int point,
                           RangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

      typedef typename LocalBaseFunctionType :: BaseFunctionSetType LocalBaseFunctionSetType;

      LocalBaseFunctionType localBaseFunction = baseFunctionList_[ baseFunction ].localFunction( entity_ );
      LocalBaseFunctionSetType &localBaseFunctionSet = localBaseFunction.baseFunctionSet();
      const int numLocalBaseFunctions = localBaseFunctionSet.numBaseFunctions();
      
      phi = 0;
      for( int i = 0; i < numLocalBaseFunctions; ++i )
      {
        RangeType psi;
        localBaseFunctionSet.evaluate( i, diffVariable, quadrature, point, psi );
        phi.axpy( localBaseFunction[ i ], psi );
      }
    }

    inline int numBaseFunctions () const
    {
      assert( baseFunctionList_ != NULL );
      return baseFunctionList_->size();
    }
  };

}

#endif
