#ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
#define DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH

#include <dune/fem/space/common/dofmapper.hh>

namespace Dune
{

  template< class GridPartImp, class BaseFunctionListImp >
  class ReducedBasisMapper;



  template< class GridPartImp, class BaseFunctionListImp >
  struct ReducedBasisMapperTraits
  {
    typedef GridPartImp GridPartType;
    
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;
    
    typedef BaseFunctionListImp BaseFunctionListType;

    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType >
      DofMapperType;
    
    typedef DefaultDofMapIterator< EntityType, DofMapperType >
      DofMapIteratorType;
  };



  /** \class ReducedBasisMapper
   *  \brief provides the mapper for the reduced basis space 
   *
   *  This mapper just performs the identity mapping.
   */
  template< class GridPartImp, class BaseFunctionListImp >
  class ReducedBasisMapper
  : public DofMapperDefault
    < ReducedBasisMapperTraits< GridPartImp, BaseFunctionListImp > >
  {
  public:
    typedef ReducedBasisMapperTraits< GridPartImp, BaseFunctionListImp >
      Traits;
   
    typedef typename Traits :: BaseFunctionListType BaseFunctionListType;
    
    typedef typename Traits :: GridPartType GridPartType;
    typedef typename Traits :: EntityType EntityType;
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

  private:
    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  protected:
    const BaseFunctionListType &baseFunctionList_;

  public:
    inline explicit ReducedBasisMapper ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( baseFunctionList )
    {
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    inline DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    inline DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal(const EntityType &entity,const int localDof) const */
    inline int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      return localDof;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal(const Entity &entity,const int localDof) const */
    template< class Entity >
    inline int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      DUNE_THROW( NotImplemented, "ReducedBasisSpace cannot map entity DoFs." );
      return 0;
    }

    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    inline int maxNumDofs () const
    {
      return size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    inline int numEntityDofs ( const Entity &entity ) const
    {
      DUNE_THROW( NotImplemented, "ReducedBasisSpace cannot map entity DoFs." );
      return 0;
    }
   
    /** \copydoc Dune::DofMapper::consecutive() const */
    bool consecutive () const
    {
      return false;
    }

    /** \copydoc Dune::DofMapper::newIndex(const int hole,const int block) const */
    int newIndex ( const int hole, const int block ) const
    {
      return -1;
    }

    /** \copydoc Dune::DofMapper::newSize() const */
    int newSize () const
    {
      return size();
    }

    /** \copydoc Dune::DofMapper::numberOfHoles(const int block) const */
    int numberOfHoles ( const int block ) const
    {
      return 0;
    }

    /** \copydoc Dune::DofMapper::oldIndex(const int hole,const int block) const */
    int oldIndex ( const int hole, const int block ) const
    {
      return -1;
    }

    /** \copydoc Dune::DofMapper::size() const */
    int size () const
    {
      return baseFunctionList_.size();
    }
  };
  
};

#endif
