#ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
#define DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>

namespace Dune
{

  template< class GridPartImp, class BaseFunctionListImp >
  class ReducedBasisMapper;



  template< class GridPartImp, class BaseFunctionListImp >
  struct ReducedBasisMapperTraits
  {
    typedef GridPartImp GridPartType;
    
    typedef typename GridPartType::template Codim< 0 >::EntityType ElementType;
    
    typedef BaseFunctionListImp BaseFunctionListType;

    typedef ReducedBasisMapper< GridPartType, BaseFunctionListType >
      DofMapperType;
    
    typedef DefaultDofMapIterator< ElementType, DofMapperType > DofMapIteratorType;
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
    explicit ReducedBasisMapper ( const BaseFunctionListType &baseFunctionList )
    : baseFunctionList_( baseFunctionList )
    {}

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType
        ( DofMapIteratorType :: endIterator, entity, *this );
    }

    /** \copydoc Dune::DofMapper::mapEach */
    template< class Functor >
    void mapEach ( const ElementType &element, Functor f ) const
    {
      const int n = size();
      for( int i = 0; i < n; ++i )
        f( i, i );
    }

    /** \copydoc Dune::DofMapper::mapToGlobal(const EntityType &entity,const int localDof) const */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      return localDof;
    }

    /** \copydoc Dune::DofMapper::mapEntityDofToGlobal(const Entity &entity,const int localDof) const */
    template< class Entity >
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      DUNE_THROW( NotImplemented, "ReducedBasisSpace cannot map entity DoFs." );
      return 0;
    }

    /** \copydoc Dune::DofMapper::numDofs(const ElementType &element) const */
    int numDofs ( const ElementType &element ) const
    {
      return size();
    }

    /** \copydoc Dune::DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return size();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const Entity &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
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
  
} // namespace Dune

#endif // #ifndef DUNE_FEM_REDUCEDBASISSPACE_MAPPER_HH
