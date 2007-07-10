#ifndef DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH
#define DUNE_FEM_REDUCEDBASISSPACE_REDUCEDBASISSPACE_HH

#include "basefunctionset.hh"
#include "mapper.hh"

namespace Dune
{

  template< class BaseFunctionImp >
  class ReducedBasisSpace
  {
  public:
    typedef BaseFunctionImp BaseFunctionType;

  private:
    typedef ReducedBasisSpace< BaseFunctionType > ThisType;

  public:
    typedef typename BaseFunctionType :: FunctionSpaceType
      BaseFunctionSpaceType;

    typedef typename BaseFunctionSpaceType :: GridPartType GridPartType;
    typedef typename BaseFunctionSpaceType :: IndexSetType IndexSetType;
    typedef typename BaseFunctionSpaceType :: IteratorType IteratorType;

    typedef ReducedBasisBaseFunctionSet< BaseFunctionType >
      BaseFunctionSetType;

    typedef typename BaseFunctionSetType :: BaseFunctionListType BaseFunctionListType;

    typedef ReducedBasisMapper< BaseFunctionListType >
      MapperType;

  protected:
    const BaseFunctionSpaceType &baseFunctionSpace_;
    const BaseFunctionListType baseFunctionList_;
    const MapperType &mapper_;

  public:
    inline ReducedBasisSpace ( const BaseFunctionSpaceType &baseFunctionSpace )
    : baseFunctionSpace_( baseFunctionSpace ),
      baseFunctionList_( NULL ),
      mapper( baseFunctionList_ )
    {
    }

    inline ~ReducedBasisSpace ()
    {
      unsigned int size = baseFunctionList_.size();
      for( unsigned int i = 0; i < size; ++i )
        delete baseFunctionList_[ i ];
    }

    inline addBaseFunction ( const BaseFunctionType &baseFunction )
    {
      BaseFunctionType *f = new BaseFunctionType( baseFunction );
      baseFunctionList_.push_back( f );
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

    inline int polynomOrder () const DUNE_DEPRECATED
    {
      return order();
    }

    //! begin iterator
    inline IteratorType begin () const
    {
      return baseFunctionSpace_.begin();
    }

    //! end iterator
    inline IteratorType end () const
    {
      return baseFunctionSpace_.end();
    }

    //! provide access to the base function set for an entity
    template< class EntityType >
    inline const BaseFunctionSetType
      baseFunctionSet( const EntityType &entity ) const
    {
      return BaseFunctionSetType( baseFunctionList_, entity );
    }

    //! get dimension of value
    inline int dimensionOfValue () const
    {
      return baseFunctionSpace_.dimensionOfValue;
    }

    //! obtain the associated grid
    inline const GridType &grid () const
    {
      return BaseFunctionSpace_.grid();
    }
    
    //! obtain the associated grid partition
    inline const GridPartType &gridPart () const
    {
      return baseFunctionSpace_.gridPart();
    }

    //! obtain the associated grid partition
    inline GridPartType &gridPart ()
    {
      return baseFunctionSpace_.gridPart();
    }

    //! obtain the associated index set
    inline const IndexSetType &indexSet () const
    {
      return baseFunctionSpace_.indexSet();
    }

    //! obtain the DoF mapper of this space
    inline MapperType &mapper () const
    {
      assert( mapper_ != NULL );
      return *mapper_;
    }

    //! map local DoF number to global DoF number
    template< class EntityType >
    inline int mapToGlobal( EntityType &entity, int localDof ) const
    {
      return mapper_->mapToGlobal( entity, localDof );
    }

    //! number of DoFs in the function space
    inline int size () const
    {
      return mapper_->size();
    }
  };
  
}

#endif
