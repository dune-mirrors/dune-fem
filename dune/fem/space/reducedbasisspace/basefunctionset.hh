#ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
#define DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/space/basefunctions/basefunctionsetinterface.hh>

namespace Dune
{

  template< class BaseFunctionImp >
  class ReducedBasisBaseFunctionSet;


  
  template< class BaseFunction >
  class ReducedBasisBaseFunctionSetTraits
  {
    typedef ReducedBasisBaseFunctionSetTraits< BaseFunction > ThisType;

  public:
    typedef ReducedBasisBaseFunctionSet< BaseFunction > BaseFunctionSetType;

    typedef BaseFunction BaseFunctionType;
    
    typedef typename BaseFunctionType::DiscreteFunctionSpaceType BaseFunctionSpaceType;
     
    typedef typename BaseFunctionSpaceType::BaseFunctionSetType::FunctionSpaceType
      FunctionSpaceType;
  };



  /** \class ReducedBasisBaseFunctionSet 
   *  \brief The ReducedBasisBaseFunctionSet class provides  
   *
   *  This class is needed to build the space and provides the functionality of
   *  the space for example the jacobian method is implemented here 
   */
  template< class BaseFunction >
  class ReducedBasisBaseFunctionSet
  : public BaseFunctionSetDefault< ReducedBasisBaseFunctionSetTraits< BaseFunction > >
  {
    typedef ReducedBasisBaseFunctionSet< BaseFunction > ThisType;
    typedef BaseFunctionSetDefault< ReducedBasisBaseFunctionSetTraits< BaseFunction> >
      BaseType;

  public:
    typedef ThisType BaseFunctionSetType;

    typedef ReducedBasisBaseFunctionSetTraits< BaseFunction > Traits;

    typedef typename Traits::BaseFunctionType BaseFunctionType;
    typedef typename Traits::BaseFunctionSpaceType BaseFunctionSpaceType;

    typedef typename BaseFunctionType::LocalFunctionType LocalBaseFunctionType;

    typedef typename BaseFunctionSpaceType::GridPartType GridPartType;
   
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;
    
    typedef typename FunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename FunctionSpaceType::RangeFieldType RangeFieldType;

    typedef typename FunctionSpaceType::DomainType DomainType;
    typedef typename FunctionSpaceType::RangeType RangeType;

    typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

    static const int dimDomain = FunctionSpaceType::dimDomain;
    static const int dimRange = FunctionSpaceType::dimRange;

    typedef DynamicArray< BaseFunctionType* > BaseFunctionListType;

  private:
    typedef typename GridPartType::GridType::template Codim< 0 >::Entity ElementType;

  public:
    using BaseType::evaluate;

    /** \brief default constructor */
    ReducedBasisBaseFunctionSet ()
    : baseFunctionList_( 0 ),
      entity_( 0 )
    {}
   
    /** constructor
     *
     *  This constructor initializes the base function set, but does not bind
     *  it to an entity.
     *
     *  \param[in]  baseFunctionList  array containing the discrete functions
     *                                to be used as base functions
     */
    explicit ReducedBasisBaseFunctionSet ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( &baseFunctionList ),
      entity_( 0 )
    {}

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
    ReducedBasisBaseFunctionSet ( const BaseFunctionListType &baseFunctionList,
                                  const ElementType &entity )
    : baseFunctionList_( &baseFunctionList ),
      entity_( &entity )
    {}
    
    /** \copydoc Dune::BaseFunctionSetInterface::evaluate(const int baseFunction,const FieldVector<int,diffOrd> &diffVariable,const PointType &x,RangeType &phi) const */
    template< int diffOrd, class PointType >
    void evaluate ( const int baseFunction,
                    const FieldVector< int, diffOrd > &diffVariable,
                    const PointType &x,
                    RangeType &phi ) const;

    /** \copydoc Dune::BaseFunctionSetInterface::evaluateSingle(const int baseFunction,const PointType &x,const RangeType &psi) const */
    template< class PointType >
    RangeFieldType evaluateSingle ( const int baseFunction,
                                    const PointType &x,
                                    RangeType &psi ) const;
   
    /** \brief obtain the entity, this base function set belongs to */
    const ElementType &entity () const
    {
      assert( entity_ != 0 );
      return *entity_;
    }

    GeometryType geometryType () const
    {
      return entity().geometry().type();
    }

    /** \copydoc Dune::BaseFunctionSetInterface::numBaseFunctions */
    int numBaseFunctions () const
    {
      assert( baseFunctionList_ != 0 );
      return baseFunctionList_->size();
    }

  protected:
    const BaseFunctionListType *baseFunctionList_;
    const ElementType *entity_;
  };

}

#include "basefunctionset_inline.hh"

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_BASEFUNCTIONSET_HH
