#ifndef DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH

#include "basefunctionset.hh"
#include "mapper.hh"

#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

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
 *  have an underlying arbitrary space. Consequently they inhert most of the 
 *  structure from this space.
 *  Initially the space is empty and by using the add function you can bulid this space 
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
  //! constructor the underlying basis is the argument
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
      clear();
    }

    /** \brief add a base function to the reduced basis space
     *
     *  \note After adding the base function, you can safely remove or
     *        overwrite it. The ReducedBasisSpace make itself a copy.
     * 
     *  \param[in]  baseFunction  base function to add to the reduced basis
     *                            space
     */
    inline void addBaseFunction ( const BaseFunctionType &baseFunction )
    {
      BaseFunctionType *f = new BaseFunctionType( baseFunction );
      assert( f != NULL );
      //baseFunctionList_.push_back( f );
      baseFunctionList_.append( f );
    }

    /** \brief access a base function within the reduced basis space
     *
     *  \param[in]  i  number of the base function to access
     * 
     *  \returns a constant reference to the i-th base function
     */
    inline const BaseFunctionType &baseFunction ( unsigned int i ) const
    {
      return *(baseFunctionList_[ i ]);
    }

    inline const BaseFunctionSpaceType &baseFunctionSpace () const
    {
      return baseFunctionSpace_;
    }

    /** \brief remove all base functions from the reduced basis space
     */
    inline void clear ()
    {
      unsigned int size = baseFunctionList_.size();
      for( unsigned int i = 0; i < size; ++i )
        delete baseFunctionList_[ i ];
      baseFunctionList_.resize( 0 );
    }
    
    /** \brief obtain number of base functions within the reduced basis space
     */
    inline unsigned int numBaseFunctions () const
    {
      return baseFunctionList_.size();
    }

    //! are the functions continuous?
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
    inline bool multipleBaseFunctionSets () const
    {
      return true;
    }
  };
  
}

#include "reducedbasisspace_inline.hh"

#endif
