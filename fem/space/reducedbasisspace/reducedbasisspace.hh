#ifndef DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH

#include "basefunctionset.hh"
#include "mapper.hh"

namespace Dune
{

  template< class BaseFunctionImp >
  class ReducedBasisSpace;



/*======================================================================*/
/*!
 *  \class ReducedBasisSpaceTraits 
 *  \brief The ReducedBasisSpaceTraits class provides  the traits for the RBspace
 *
 *  many typedefs
 *
 */
/*======================================================================*/
  template< class BaseFunctionImp >
  class ReducedBasisSpaceTraits
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

  private:
    typedef ReducedBasisSpaceTraits< BaseFunctionType > ThisType;

  public:
    typedef ReducedBasisSpace< BaseFunctionType > DiscreteFunctionSpaceType;

    typedef typename BaseFunctionType :: FunctionSpaceType BaseFunctionSpaceType;

    typedef typename BaseFunctionSpaceType :: DomainType DomainType;
    typedef typename BaseFunctionSpaceType :: RangeType RangeType;
    typedef typename BaseFunctionSpaceType :: JacobianRangeType JacobianRangeType;
    
    typedef typename BaseFunctionSpaceType :: DomainFieldType DomainFieldType;
    typedef typename BaseFunctionSpaceType :: RangeFieldType RangeFieldType;

    typedef typename BaseFunctionSpaceType :: FunctionSpaceType FunctionSpaceType;
    typedef typename BaseFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseFunctionSpaceType :: GridType GridType;
    typedef typename BaseFunctionSpaceType :: IndexSetType IndexSetType;
    typedef typename BaseFunctionSpaceType :: IteratorType IteratorType;

    typedef ReducedBasisBaseFunctionSet< BaseFunctionType >
      BaseFunctionSetType;

    typedef typename BaseFunctionSetType :: BaseFunctionListType BaseFunctionListType;

    typedef ReducedBasisMapper< BaseFunctionListType >
      MapperType;
  };


  

/*======================================================================*/
/*!
 *  \class ReducedBasisSpace 
 *  \brief The ReducedBasisSpace class provides the space for RB simulations 
 *
 *  The basis consists of discrete functions as basis functions. These discrete functions 
 *  have an underlying space, the Lagrange space. Consequently they inhert most of the 
 *  structure from the Lagrange space.
 *  Initially the space is empty and be using the add function you can bulid this space 
 *  and discrete functions. 
 *  
 *
 */
/*======================================================================*/
  template< class BaseFunctionImp >
  class ReducedBasisSpace
  : public DiscreteFunctionSpaceDefault< ReducedBasisSpaceTraits< BaseFunctionImp > >
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

    typedef ReducedBasisSpaceTraits< BaseFunctionType > Traits;

  private:
    typedef ReducedBasisSpace< BaseFunctionType > ThisType;
    typedef DiscreteFunctionSpaceDefault< Traits > BaseType;

  public:
    typedef typename Traits :: BaseFunctionSpaceType BaseFunctionSpaceType;
    typedef typename Traits :: FunctionSpaceType FunctionSpaceType;
    typedef typename Traits :: GridPartType GridPartType;
    typedef typename Traits :: GridType GridType;
    typedef typename Traits :: IndexSetType IndexSetType;
    typedef typename Traits :: IteratorType IteratorType;

    typedef typename Traits :: BaseFunctionSetType BaseFunctionSetType;
    typedef typename Traits :: BaseFunctionListType BaseFunctionListType;
    typedef typename Traits :: MapperType MapperType;

    enum { polynomialOrder = BaseFunctionSpaceType :: polynomialOrder };

  protected:
    BaseFunctionSpaceType &baseFunctionSpace_;
    BaseFunctionListType baseFunctionList_;
    mutable MapperType mapper_;

  public:
  //! constructor the underlying lagrange basis is the argument
    inline explicit ReducedBasisSpace ( BaseFunctionSpaceType &baseFunctionSpace )
    : BaseType( baseFunctionSpace.gridPart() ),
      baseFunctionSpace_( baseFunctionSpace ),
      baseFunctionList_(),
      mapper_( baseFunctionList_ )
    {
    }

  //! destructor to release the pointer of each entry
    inline ~ReducedBasisSpace ()
    {
      unsigned int size = baseFunctionList_.size();
      for( unsigned int i = 0; i < size; ++i )
        delete baseFunctionList_[ i ];
    }

   //! this method is used to create the reduced basis space by adding discrete functions
    inline void addBaseFunction ( const BaseFunctionType &baseFunction )
    {
      BaseFunctionType *f = new BaseFunctionType( baseFunction );
      assert( f != NULL );
      //baseFunctionList_.push_back( f );
      baseFunctionList_.append( f );
    }
    
    //! are the function continuous?
    inline bool continuous () const
    {
      return baseFunctionSpace_.continuous();
    }
    
    //! get the polynomial order of this discrete function space
    inline int order () const
    {
      return baseFunctionSpace_.order();
    }

    //! provide access to the base function set for an entity
    template< class EntityType >
    inline const BaseFunctionSetType baseFunctionSet( const EntityType &entity ) const
    {
      return BaseFunctionSetType( baseFunctionList_, entity );
    }

    //! get dimension of value
    inline int dimensionOfValue () const
    {
      return baseFunctionSpace_.dimensionOfValue;
    }

    //! obtain the DoF mapper of this space
    inline MapperType &mapper () const
    {
      return mapper_;
    }

    //! are there multiple base function sets per geometry type?
    inline bool multipleBaseFunctionSets() const
    {
      return true;
    }

  };
  
}

#endif
