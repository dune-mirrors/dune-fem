#ifndef DUNE_FEM_CODIMENSIONMAPPER_HH
#define DUNE_FEM_CODIMENSIONMAPPER_HH

#include <dune/fem/space/common/dofmapper.hh>

namespace Dune
{

  template< class GridPartImp, int codim >
  class CodimensionMapper;

  template <class EntityType, class DofMapperType> 
  class CodimensionDofMapIterator;


  template< class GridPartImp, int codim >
  struct CodimensionMapperTraits
  {
    typedef GridPartImp GridPartType;

    // we still need entities of codimension 0 here 
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    
    typedef CodimensionMapper< GridPartType, codim >  DofMapperType;

    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;
  };


  template< class GridPartImp >
  struct CodimensionMapperTraits< GridPartImp, 0 >
  {
    typedef GridPartImp GridPartType;

    // we still need entities of codimension 0 here 
    typedef typename GridPartType :: template Codim< 0 > :: IteratorType :: Entity
      EntityType;

    typedef typename GridPartType :: IndexSetType IndexSetType;
    
    typedef CodimensionMapper< GridPartType, 0 >  DofMapperType;

    typedef CodimensionDofMapIterator<EntityType, DofMapperType>  DofMapIteratorType;
  };



  //***************************************************************************
  //
  //!  The CodimensionMapper for mapping of local dof numbers to global dof numbers, 
  //!  i.e. the entry in the vector of unknowns
  //
  //***************************************************************************
  template< class GridPartImp, int cdim >
  class CodimensionMapper
  : public DofMapperDefault< CodimensionMapperTraits< GridPartImp, cdim > >
  {
  public:
    //! codimension that is mapped 
    enum { codimension = cdim };

    typedef CodimensionMapperTraits< GridPartImp, codimension> Traits;

    typedef typename Traits :: EntityType EntityType;

    typedef typename Traits :: GridPartType GridPartType;

    typedef typename Traits :: IndexSetType IndexSetType;

    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

  private:
    typedef CodimensionMapper< GridPartType, codimension > ThisType;
    typedef DofMapperDefault< Traits > BaseType;

  protected:
    // index set of grid, i.e. the element indices 
    const IndexSetType &indexSet_;

    // max number of local dofs 
    const int maxNumberOfDofs_;
  public:
    //! Constructor 
    inline CodimensionMapper( const GridPartType &gridPart,
                              const int maxDofs )
    : indexSet_( gridPart.indexSet() ),
      maxNumberOfDofs_( maxDofs )
    {}

    /** \copydoc DofMapper::contains */
    bool contains ( const int codim ) const 
    {
      return ( codimension == codim );
    }

    //! return size of function space 
    //! see dofmanager.hh for definition of IndexSet, which 
    //! is a wrapper for en.index 
    /** \copydoc DofMapper::size */
    int size () const
    {
      // return number of dofs for codimension 
      return indexSet_.size( codimension );
    }

    /** \copydoc Dune::DofMapper::begin(const EntityType &entity) const */
    DofMapIteratorType begin ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    /** \copydoc Dune::DofMapper::end(const EntityType &entity) const */
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    /** \copydoc DofMapper::mapToGlobal */
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      // we only have one local dof 
      assert( localDof < maxNumDofs() );
      if ( codimension == 0 ) 
        return indexSet_.index( entity );
      else 
        return indexSet_.subIndex( entity, localDof, codimension );
    }

    /** \copydoc DofMapper::maxNumDofs() const */
    int maxNumDofs () const
    {
      return maxNumberOfDofs_;
    }

    using BaseType::numDofs;

    /** \copydoc Dune::DofMapper::numDofs(const EntityType &entity) const */
    int numDofs ( const EntityType &entity ) const
    {
      return entity.template count< codimension >();
    }

    /** \copydoc Dune::DofMapper::numEntityDofs(const EntityType &entity) const */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      // check for codimension equality 
      return ( ( Entity::codimension - codimension ) == 0 ) ? 1 : 0;
    }
   
    /** \copydoc DofMapper::oldIndex */
    int oldIndex (const int hole, int ) const
    {
      // forward to index set 
      return indexSet_.oldIndex( hole, codimension ) ;
    }

    /** \copydoc DofMapper::newIndex */
    int newIndex (const int hole, int ) const
    {
      // forward to index set 
      return indexSet_.newIndex( hole, codimension );
    }

    /** \copydoc DofMapper::numberOfHoles */
    int numberOfHoles ( const int ) const
    {
      // this index set works only for codim = 0 at the moment
      return indexSet_.numberOfHoles( codimension );
    }

    /** \copydoc DofMapper::consecutive */
    bool consecutive() const 
    {
      return BaseType::checkConsecutive(indexSet_);
    }
  };


  template <class EntityType, class DofMapperType> 
  class CodimensionDofMapIterator
  {
  public:
    enum IteratorType { beginIterator, endIterator };

  private:
    typedef CodimensionDofMapIterator<EntityType, DofMapperType>   ThisType;
    
  protected:
    const int baseIndex_;
    int dof_;
    
  public:
    inline CodimensionDofMapIterator ( const IteratorType type,
                                       const EntityType& entity,
                                       const DofMapperType& mapper)
    : baseIndex_( (type == endIterator) ? -1 : 
                    mapper.mapToGlobal( entity, 0 ) ),
      dof_( (type == endIterator) ? mapper.numDofs( entity ) : 0 )
    {}

    inline CodimensionDofMapIterator ( const ThisType &other )
    : baseIndex_( other.baseIndex_ ),
      dof_( other.dof_ )
    {}

    inline ThisType &operator++ ()
    {
      ++dof_;
      return *this;
    }

    inline bool operator== ( const ThisType &other ) const
    {
      return dof_ == other.dof_;
    }

    inline bool operator!= ( const ThisType &other ) const
    {
      return dof_ != other.dof_;
    }

    inline int local () const
    {
      return dof_;
    }

    inline int global () const
    {
      return baseIndex_ + dof_;
    }
  };

} // end namespace Dune

#endif
