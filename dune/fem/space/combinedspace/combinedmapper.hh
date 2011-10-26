#ifndef DUNE_COMBINEDMAPPER_HH
#define DUNE_COMBINEDMAPPER_HH

#include <dune/fem/space/mapper/dofmapper.hh>
#include <dune/fem/space/common/dofmanager.hh>



namespace Dune
{

  //! forward declaration
  template< class Grid,  class Mapper1, class Mapper2 >
  class CombinedMapper;


  //! CombinedDofMapIterator 
  template<class EntityType, class DofMapperType >
  class CombinedDofMapIterator;

  //! Traits 
  template< class Grid, class Mapper1, class Mapper2 >
  struct CombinedMapperTraits
  {

    // we still need entities of codimension 0 here 
    typedef typename Mapper1 :: EntityType EntityType;

    static const int polynomialOrder1 = Mapper1 :: polynomialOrder;
    static const int polynomialOrder2 = Mapper2 :: polynomialOrder;
    static const int polynomialOrder = ( polynomialOrder1 > polynomialOrder2 ) ? polynomialOrder1 : polynomialOrder2;

//    typedef typename Mapper1 :: IndexSetType IndexSetType;

    typedef CombinedMapper< Grid, Mapper1, Mapper2 >  DofMapperType;

    // Later on write an new dofIterator
    typedef DefaultDofMapIterator< EntityType, DofMapperType > DofMapIteratorType;
  };



  template< class Grid, class Mapper1, class Mapper2 >
  class CombinedMapper
  : public DofMapperDefault< CombinedMapperTraits< Grid, Mapper1, Mapper2 > >
  {
    typedef CombinedMapper< Grid, Mapper1, Mapper2 > ThisType;
    typedef DofMapperDefault< CombinedMapperTraits< Grid, Mapper1, Mapper2 > > BaseType;
    typedef Grid GridType;

    typedef Mapper1 MapperType1;
    typedef Mapper2 MapperType2;

    public:
    typedef typename BaseType :: Traits Traits;
    typedef typename BaseType :: EntityType EntityType;

//    typedef typename Traits :: GridPartType GridPartType;
//    typedef typename Traits :: IndexSetType IndexSetType;
    typedef typename Traits :: DofMapIteratorType DofMapIteratorType;

    //! type of the underlying grid
//    typedef typename GridPartType::GridType GridType;
    //! type of coordinates within the grid
//    typedef typename GridType::ctype FieldType;
    //! dimension of the grid
//    static const int dimension = GridType::dimension;
    //! order of the Lagrange polynoms
    static const int polynomialOrder1 = Traits :: polynomialOrder1;
    static const int polynomialOrder2 = Traits :: polynomialOrder2;
    static const int polynomialOrder = Traits :: polynomialOrder;

    //! type of the DoF manager
    typedef DofManager< GridType > DofManagerType;

    public:
    //! constructor
    CombinedMapper( const GridType &grid,  const MapperType1 &mapper1, const MapperType2& mapper2 )
      : // BaseType( gridPart ), // hmmm
        dm_( DofManagerType :: instance( grid ) ), 
        mapper1_( mapper1 ),
        mapper2_( mapper2 ),
        globalOffset_( mapper1.size() ),
        oldGlobalOffset_( -1 )
    {
      dm_.addIndexSet( *this );
    }

    ~CombinedMapper() 
    {
      dm_.removeIndexSet( *this );
    }
    
    //! copy constructor
    CombinedMapper( const ThisType &other )
      : dm_( other.dm_ ),
        mapper1_( other.mapper1_ ),
        mapper2_( other.mapper2_ ),
        globalOffset_( other.globalOffset_ ),
        oldGlobalOffset_( other.oldGlobalOffset_ )
    {
      dm_.addIndexSet( *this );
    }


    int size () const 
    {
      return mapper1_.size() + mapper2_.size();
    }

    // needs to be reimpl
    DofMapIteratorType begin ( const EntityType &entity ) const   
    {
      return DofMapIteratorType( DofMapIteratorType::beginIterator, entity, *this );
    }
    
    DofMapIteratorType end ( const EntityType &entity ) const
    {
      return DofMapIteratorType( DofMapIteratorType::endIterator, entity, *this );
    }

    bool contains ( const int codim ) const    
    {
      return ( mapper1_.contains( codim ) || mapper2_.contains( codim ) );
    }
    
    int mapToGlobal ( const EntityType &entity, const int localDof ) const
    {
      assert( mapper1_.size() == globalOffset_ );
      const int localOffset = mapper1_.numDofs(entity);

      int index;

      if( localDof - localOffset  < 0 )
        index = mapper1_.mapToGlobal( entity, localDof );
      else
        index = mapper2_.mapToGlobal( entity, localDof - localOffset ) + globalOffset_;
        
      assert( (0 <= index) && (index < size()) );
      return index;
    }
    
    template< class Entity > 
    int mapEntityDofToGlobal ( const Entity &entity, const int localDof ) const
    {
      assert( mapper1_.size() == globalOffset_ );
      const int numEntityDofs =  mapper1_.numEntityDofs(entity);

      if( localDof - numEntityDofs < 0 )
        return mapper1_.mapEntityDofToGlobal( entity, localDof );
      else
        return mapper2_.mapEntityDofToGlobal( entity, localDof - numEntityDofs ) + globalOffset_;
    }
    
    /** \brief obtain maximal number of DoFs on one entity
     */
    int maxNumDofs () const
    {
      return mapper1_.maxNumDofs() + mapper2_.maxNumDofs();
    }

    /** \brief obtain number of DoFs on an entity
     * 
     *  \param[in]  entity  entity of codimension 0
     *  
     *  \returns number of DoFs on the entity
     */
    int numDofs ( const EntityType &entity ) const
    {
      int nDofs = mapper1_.numDofs( entity ) + mapper2_.numDofs( entity );
//      assert( (nDofs >=0) && (nDofs < size()) );
      return nDofs;
    }

    /** \brief obtain number of DoFs actually belonging to an entity
     *
     *  In contrast to numDofs, this method returns the number of DoFs actually
     *  associated with an entity (usually a subentity). We have the following
     *  relation for an entity \f$E\f$ of codimension 0:
     *  \f[
     *  \mathrm{numDofs}( E ) = \sum_{e \subset E} \mathrm{numEntityDofs}( e ),
     *  \f]
     *  where \f$\subset\f$ denotes the subentity relation.
     * 
     *  \param[in]  entity  entity of codimension
     *  
     *  \returns number of DoFs on the entity
     */
    template< class Entity >
    int numEntityDofs ( const Entity &entity ) const
    {
      return mapper1_.numEntityDofs( entity ) + mapper2_.numEntityDofs( entity );
    }

    /** \brief Check, whether the data in a codimension has fixed size */
    bool fixedDataSize ( int codim ) const
    {
      return mapper1_.fixedDataSize( codim ) && mapper2_.fixedDataSize( codim );
    }

    /** \brief return number of holes for data block */
    int numberOfHoles(const int block) const 
    {
      const int numBlock1 = mapper1_.numBlocks(); 
      if ( block - numBlock1 < 0 ) 
        return mapper1_.numberOfHoles( block );
      else
        return mapper2_.numberOfHoles( block - numBlock1 );
    }
    
    /** \brief return old index of hole for data block (with resprect to new offset) */
    int oldIndex (const int hole, const int block) const 
    { 
      const int numBlock1 = mapper1_.numBlocks(); 
      if ( block - numBlock1 < 0 ) 
        return mapper1_.oldIndex( hole, block );
      else
        return mapper2_.oldIndex( hole, block - numBlock1 );
    }
      
    /** \brief return new index of hole for data block (with resprect to new offset) */
    int newIndex (const int hole, const int block) const 
    { 
      const int numBlock1 = mapper1_.numBlocks(); 
      if ( block - numBlock1 < 0 ) 
        return mapper1_.newIndex( hole, block );
      else
        return mapper2_.newIndex( hole, block - numBlock1 );
    }

    /** \brief return true if compress will affect data */
    bool consecutive () const 
    {
      return mapper1_.consecutive() && mapper2_.consecutive();
    }

    /** \brief return old offsets for given block */
    int oldOffSet(const int block) const
    {
      const int numBlock1 = mapper1_.numBlocks();
      if ( block - numBlock1 < 0 )
        return mapper1_.oldOffSet( block );
      else
        return mapper2_.oldOffSet( block - numBlock1 ) + oldGlobalOffset_;
    }

    /** \brief return current offsets for given block */
    int offSet(const int block) const
    {
      assert( globalOffset_ == mapper1_.size() );
      const int numBlock1 = mapper1_.numBlocks(); 
      if ( block - numBlock1 < 0 ) 
        return mapper1_.offSet( block );
      else
        return mapper2_.offSet( block - numBlock1 ) + globalOffset_;
    }

    /** \brief return number of supported blocks */
    int numBlocks() const
    {
      return mapper1_.numBlocks() + mapper2_.numBlocks();
    }

    /** \brief return number of supported blocks */
    void resize ()
    {
      oldGlobalOffset_ = globalOffset_;
      globalOffset_ = mapper1_.size();      
    }

    /** \brief return number of supported blocks */
    bool compress ()
    {
      resize();
      return true;
    }


    void read_xdr ( const char *filename, int timestep )
    {
      resize();
    }

    void write_xdr ( const char *filename, int timestep )
    {}

    void insertEntity ( const EntityType &entity )
    {
      resize();
    }

    void removeEntity ( const EntityType &entity )
    {
      resize();
    }


  protected:
    DofManagerType &dm_;
    const MapperType1 &mapper1_;
    const MapperType2 &mapper2_;

    int globalOffset_;
    int oldGlobalOffset_;
  };


  template<class EntityType, class DofMapperType >
  class CombinedDofMapIterator
  {
    typedef CombinedDofMapIterator< EntityType, DofMapperType > ThisType;

  public:
    enum IteratorType { beginIterator, endIterator };

    CombinedDofMapIterator ( const IteratorType type,
                             const EntityType& entity,
                             const DofMapperType& mapper)
    : baseIndex_( (type == endIterator) ? -1 : 
                   mapper.mapToGlobal( entity, 0 ) ),
      dof_( (type == endIterator) ? mapper.numDofs( entity ) : 0 )
    {}

    CombinedDofMapIterator ( const ThisType &other )
    : baseIndex_( other.baseIndex_ ),
      dof_( other.dof_ )
    {}

    const ThisType &operator++ ()
    {
      ++dof_;
      return *this;
    }

    bool operator== ( const ThisType &other ) const
    {
      return dof_ == other.dof_;
    }

    bool operator!= ( const ThisType &other ) const
    {
      return dof_ != other.dof_;
    }

    int local () const
    {
      return dof_;
    }

    int global () const
    {
      return baseIndex_ + dof_;
    }

  protected:
    const int baseIndex_;
    int dof_;

  };

#if 0
  // ComobinedMapperSingletonFactory
  // ---------------------------------

  template<class Grid, class Mapper1, class Mapper2 >
  struct CombinedMapperSingletonFactory
  {
    typedef CombinedMapper< Grid, Mapper1, Mapper2 > Object;

    struct Key
    {
      Key ( const Mapper1 &mp1, const Mapper2 &mp2 )
      : mp1_( mp1 ),
        mp2_( mp2 )
      {}

      bool operator== ( const Key &other )
      {
        return ((&mp1_ == &other.mp1_) && (&mp2_ == &other.mp2_))  ;
      }

      const Mapper1 &mp1_;
      const Mapper2 &mp2_;
    };

    //! create new mapper  
    static Object *createObject ( const Key &key )
    {
      return new Object ( key.mp1_, key.mp2_ );
    }

    //! delete mapper object 
    static void deleteObject ( Object *object )
    {
      delete object;
    }
  };

#endif
}
#endif //header guards
