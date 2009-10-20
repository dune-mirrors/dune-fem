#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH

#include <dune/fem/storage/array.hh>
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



  /** \class ReducedBasisBaseFunctionSet 
   *  \brief The ReducedBasisBaseFunctionSet class provides  
   *
   *  This class is needed to build the space and provides the functionality of
   *  the space for example the jacobian method is implemented here 
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
    using BaseType :: evaluateSingle;

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

    enum { dimDomain = FunctionSpaceType :: dimDomain };
    enum { dimRange = FunctionSpaceType :: dimRange };

    typedef DynamicArray< BaseFunctionType* > BaseFunctionListType;

  private:
    typedef typename GridPartType :: GridType :: template Codim< 0 > :: Entity
      EntityCodim0Type;

  protected:
    const BaseFunctionListType *baseFunctionList_;
    int numBaseFunctions_;
    
    const EntityCodim0Type *entity_;

  public:
    /** \brief default constructor */
    inline ReducedBasisBaseFunctionSet ()
    : baseFunctionList_( NULL ),
      entity_( NULL )
    {
    }
   
    /** constructor
     *
     *  This constructor initializes the base function set, but does not bind
     *  it to an entity.
     *
     *  \param[in]  baseFunctionList  array containing the discrete functions
     *                                to be used as base functions
     */
    inline explicit ReducedBasisBaseFunctionSet
      ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( &baseFunctionList ),
      entity_( NULL )
    {
    }

    /** constructor
     *
     *  This constructor initializes the base function set and binds it to an
     *  entity.
     *
     *  \param[in]  baseFunctionList  array containing the discrete functions
     *                                to be used as base functions
     *  \param[in]  entity            entity (of codim 0) to bind the base
     *                                function set to
     */
    inline ReducedBasisBaseFunctionSet ( const BaseFunctionListType &baseFunctionList,
                                         const EntityCodim0Type &entity )
    : baseFunctionList_( &baseFunctionList ),
      entity_( &entity )
    {
    }
    
    /** \brief copy constructor
     *
     *  \param[in]  other  base function set to copy
     */
    inline ReducedBasisBaseFunctionSet ( const ThisType &other )
    : baseFunctionList_( other.baseFunctionList_ ),
      entity_( other.entity_ )
    {
    }

    /** \brief copy another ReducedBasisBaseFunctionSet
     *
     *  \param[in]  other  base function set to copy
     */
    inline ThisType &operator= ( const ThisType &other )
    {
      baseFunctionList_ = other.baseFunctionList_;
      entity_ = other.entity_;
      return *this;
    }

    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<deriType,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    inline void evaluate ( const int baseFunction,
                           const FieldVector< deriType, diffOrd > &diffVariable,
                           const PointType &x,
                           RangeType &phi ) const;

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    inline RangeFieldType evaluateSingle ( const int baseFunction,
                                           const PointType &x,
                                           RangeType &psi ) const;
   
    /** \brief obtain the entity, this base function set belongs to */
    inline const EntityCodim0Type &entity () const
    {
      assert( entity_ != NULL );
      return *entity_;
    }

    inline GeometryType geometryType () const
    {
      return entity().geometry().type();
    }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    inline int numBaseFunctions () const
    {
      assert( baseFunctionList_ != NULL );
      return baseFunctionList_->size();
    }
  };

}

#include "basefunctionset_inline.hh"

#endif
