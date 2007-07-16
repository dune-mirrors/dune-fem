#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH

#include <dune/fem/space/common/basefunctioninterface.hh>

namespace Dune
{

  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSet;


  
  /*! \class  ReducedBasisBaseFunctionSetTraits 
   *  \brief The  ReducedBasisBaseFunctionSetTraits class provides  typedefs
   *
   *  many typedefs
   */
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



  /*! \class ReducedBasisBaseFunctionSet 
   *  \brief The ReducedBasisBaseFunctionSet class provides  
   *
   *  this class is needed to build the space and provides the functionality of the space
   *  for example the jacobian method is implemented here 
   */
  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSet
  : public BaseFunctionSetDefault< ReducedBasisBaseFunctionSetTraits< BaseFunctionImp > >
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

    typedef ReducedBasisBaseFunctionSetTraits< BaseFunctionType > TraitsType;

  private:
    typedef ReducedBasisBaseFunctionSet< BaseFunctionType > ThisType;
    typedef BaseFunctionSetDefault< TraitsType > BaseType;

    using BaseType :: evaluate;

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

    typedef std :: vector< BaseFunctionType* > BaseFunctionListType;

  private:
    typedef typename GridPartType :: GridType :: template Codim< 0 > :: Entity EntityCodim0Type;

  protected:
    const BaseFunctionListType *baseFunctionList_;
    int numBaseFunctions_;
    
    const EntityCodim0Type *entity_;

  public:
    //! constructor
    inline ReducedBasisBaseFunctionSet ()
    : baseFunctionList_( NULL ),
      entity_( NULL )
    {
    }
   
    //! constructor with an argument as TODO
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
    
    //! copy constructor
    inline ReducedBasisBaseFunctionSet ( const ThisType &other )
    : baseFunctionList_( other.baseFunctionList_ ),
      entity_( other.entity_ )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      baseFunctionList_ = other.baseFunctionList_;
      entity_ = other.entity_;
      return *this;
    }

    //! essential method to calculate the basisfunction or their derivative of order diffOrd on a given point, &x and stores the value as RangeType in &phi 
    template< int diffOrd >
    inline void evaluate ( int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const DomainType &x,
                           RangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

      typedef typename LocalBaseFunctionType :: BaseFunctionSetType LocalBaseFunctionSetType;

      LocalBaseFunctionType localBaseFunction
        = (*baseFunctionList_)[ baseFunction ]->localFunction( *entity_ );
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

    //! essential method to calculate the basisfunction or their derivative of order diffOrd on a given point, &x and stores the value as RangeType in &phi 
    //TODO wo ist der unterschied zu obigem???  
    template< int diffOrd, class QuadratureType >
    inline void evaluate ( int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const QuadratureType &quadrature,
                           int point,
                           RangeType &phi ) const
    {
      assert( (baseFunction >= 0) && (baseFunction < numBaseFunctions()) );

      typedef typename LocalBaseFunctionType :: BaseFunctionSetType LocalBaseFunctionSetType;

      LocalBaseFunctionType localBaseFunction 
        = (*baseFunctionList_)[ baseFunction ]->localFunction( *entity_ );
      LocalBaseFunctionSetType localBaseFunctionSet = localBaseFunction.baseFunctionSet();
      const int numLocalBaseFunctions = localBaseFunctionSet.numBaseFunctions();
      
      phi = 0;
      for( int i = 0; i < numLocalBaseFunctions; ++i )
      {
        RangeType psi;
        localBaseFunctionSet.evaluate( i, diffVariable, quadrature, point, psi );
        phi.axpy( localBaseFunction[ i ], psi );
      }
    }

    //! returns the number of Discretefunctions that bulid the reduced basis space
    inline int numBaseFunctions () const
    {
      assert( baseFunctionList_ != NULL );
      return baseFunctionList_->size();
    }
  };

}

#endif
