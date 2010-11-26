#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <map>
#include <vector>

#include <dune/common/misc.hh>
#include <dune/common/forloop.hh>
#include <dune/grid/common/capabilities.hh>


namespace Dune
{
  // Extenal Forward Declarations
  // ----------------------------

  template< int dim, int dimworld >
  class AlbertaGrid;

  template< int dim, int dimworld >
  class ALUConformGrid;

  template< int dim, int dimworld >
  class ALUCubeGrid;

  template< int dim, int dimworld >
  class ALUSimplexGrid;

  template< class HostGrid, class CoordFunction, class Allocator >
  class GeometryGrid;


  /** \brief A class for storing data during an adaptation cycle.
   *
   * This container allows to store data which is to remain persistent
   * even during adaptation cycles. There is a default implementation based
   * on std::map but any grid implementation can provide a specialized implementation.
   *
   * The container expects that the class Data has a default constructor. This default
   * value is returned if an data for an entity is required which has not yet been added
   * to the container. The class also provides a size method and iterators. The
   * guarantee being that the iteration is at least over all entries in the container 
   * which are not default but the iteration could also include entries with a default
   * value. The size method returns the number of elements visited by an iteration so
   * that size() is larger or equal to the entries which are not set to a default value.
   *
   * After the grid has changed update() or reserve() should to be called.
   * In contrast to the method reserve, the method update can compress the data to
   * reduce memory consumption, whereas reserve() can be implemented not to reduce any
   * memory even if the size of the grid has changed. 
   */
  template < class Grid, class Data, class Allocator=std::allocator<Data> >  
  class PersistentContainerVector;

  /** \brief An implementation for the PersistentContainer based on a container
   * satisfying the std::vector interface and using a class providing an IndexSet
   * for storing the Data*/
  template <class Grid, class Index, class Vector>  
  class PersistentContainerVector
  {
  protected:
    typedef typename Vector::value_type Data;
    typedef Grid GridType;
    const int codim_;
    const Index& index_;
    const double overEstimate_;
    Vector data_;
    
  public:  
    //! \brief entity of codimension 0 
    typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 

    //! \brief iterator type 
    typedef typename Vector :: iterator Iterator ;
    //! \brief const iterator type 
    typedef typename Vector :: const_iterator ConstIterator ;

    //! \brief constructor creating vector container filled with default value
    //         store data on entities of given codim using index to store data in vector.
    //         The overEstimate parameter can be used to allocate more memory than
    //         required to store the data.
    PersistentContainerVector( const GridType& grid, const int codim,
                               const Index& index,
                               const double overEstimate )
      : codim_( codim )
      , index_( index )
      , overEstimate_( overEstimate )
      , data_(index.size()*overEstimate)
    {
      data_.resize(index.size());
    }

    //! \brief copy constructor 
    PersistentContainerVector( const PersistentContainerVector& other ) 
      : codim_( other.codim_ )
      , index_( other.index_ )
      , overEstimate_( other.overEstimate_ )
      , data_( other.data_ )  
    {}

    //! \brief random access to entity data with correct codimension 
    template <class Entity> 
    Data& operator [] (const Entity& entity ) 
    { 
      assert( Entity :: codimension == codim_ );
      assert( index_.index( entity ) < (typename Index::IndexType) data_.size() );
      return data_[ index_.index( entity ) ];
    }

    //! \brief random access to entity data with correct codimension 
    template <class Entity> 
    const Data& operator [] (const Entity& entity ) const
    { 
      assert( Entity :: codimension == codim_ );
      assert( index_.index( entity ) < (typename Index::IndexType) data_.size() );
      return data_[ index_.index( entity ) ];
    }

    //! \brief access for sub entity data
    Data& operator () (const ElementType& element, const int subEntity ) 
    {
      assert( index_.subIndex( element, subEntity, codim_ ) < (typename Index::IndexType) data_.size() );
      return data_[ index_.subIndex( element, subEntity, codim_ ) ];
    }

    //! \brief access for sub entity data
    const Data& operator () (const ElementType& element, const int subEntity ) const 
    {
      assert( index_.subIndex( element, subEntity, codim_ ) < (typename Index::IndexType) data_.size() );
      return data_[ index_.subIndex( element, subEntity, codim_ ) ];
    }

    //! \brief const iterator begin 
    Iterator begin() 
    {
      return data_.begin();
    }

    //! \brief const iterator begin 
    ConstIterator begin() const
    {
      return data_.begin();
    }

    //! \brief iterator end 
    Iterator end() 
    {
      return data_.end();
    }

    //! \brief const iterator end 
    ConstIterator end() const
    {
      return data_.end();
    }

    //! \brief return size of allocated data 
    size_t size() const { return data_.size(); }

    //! \brief enlarge container, compress is not necessary but could be done
    void reserve( )
    {
      if( (typename Index::IndexType) index_.size( codim_ ) > (typename Index::IndexType) data_.size() ) 
        update( );
    }

    //! \brief adjust container to correct size and set all values to default
    void clear( )
    {
      const size_t newSize = index_.size( codim_ );
      data_.resize( newSize );
      data_.clear();
    }

    //! \brief adjust container to correct size including compress 
    void update( )
    { // this could be more sophisticated (although std::vector is not stupid and
      // overestimated on its own...
      const size_t newSize = index_.size( codim_ );
      if (newSize < data_.capacity())
        data_.resize(newSize);
      else
      {
        data_.reserve(newSize*overEstimate_);
        data_.resize(newSize);
      }
    }
  };

  /** \brief An implementation for the PersistentContainer based on a container
   * satisfying the std::map interface and using a class providing an IdSet
   * for storing the Data*/
  template <class Grid, class Id, class Map>
  class PersistentContainerMap
  {
    typedef PersistentContainerMap< Grid, Id, Map > ThisType;

  protected:
    typedef typename Map :: mapped_type Data;
    typedef typename Id :: IdType  IdType;
    typedef Grid GridType;
    const GridType& grid_;
    const int codim_;
    const Id& id_;
    mutable Map data_;

    typedef typename Map :: iterator iterator ;
    typedef typename Map :: const_iterator const_iterator ;

    template <class IteratorType>
    class MyIterator
    {
      IteratorType it_;
    public: 
      MyIterator(const IteratorType& it) : it_( it ) {}
      MyIterator(const MyIterator& other) : it_( other.it_ ) {}

      bool operator == (const MyIterator& other) const { return it_ == other.it_; }
      bool operator != (const MyIterator& other) const  { return it_ != other.it_; }

      MyIterator& operator ++ () 
      {
        ++it_;
        return *this;
      }
      Data& operator * () { return (*it_).second; }
      Data* operator -> () { return &((*it_).second); }
      MyIterator& operator = (const MyIterator& other) 
      {
        it_ = other.it_;
        return *this;
      }
    };

    template< int codim , bool gridHasCodim >
    struct AdaptCodimBase
    {
      static void apply ( ThisType &container, const Data& value , const int myCodim)
      {
        if( codim == myCodim )
          container.template adaptCodim< codim > ( value );
      }
    };

    template< int codim >
    struct AdaptCodimBase< codim, false >
    {
      static void apply ( ThisType &container, const Data& value , const int myCodim)
      {
      }
    };

    template< int codim >
    struct AdaptCodim
      : public AdaptCodimBase< codim, Capabilities :: hasEntity < GridType, codim > :: v >
    {
    };

  public:  
    typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 
    typedef MyIterator< iterator > Iterator;
    typedef MyIterator< const_iterator > ConstIterator;

    //! \brief constructor creating container filled with default values.
    //
    //         Container is to be used to store data on entities of given codim using id to store data in map.
    PersistentContainerMap( const GridType& grid, const int codim, const Id& id )
      : grid_( grid )
      , codim_( codim )
      , id_( id )
      , data_()
    {
    }

    //! \brief copy constructor 
    PersistentContainerMap( const PersistentContainerMap& other ) 
      : grid_( other.grid_ )
      , codim_( other.codim_ )
      , id_( other.id_ )
      , data_( other.data_ )  
    {}

    //! \brief random access entity with correct codimension 
    template <class Entity> 
    Data& operator [] (const Entity& entity ) 
    { 
      assert( Entity :: codimension == codim_ );
      return data_[ id_.id( entity ) ];
    }

    //! \brief random access entity with correct codimension 
    template <class Entity> 
    const Data& operator [] (const Entity& entity ) const
    { 
      assert( Entity :: codimension == codim_ );
      return data_[ id_.id( entity ) ];
    }

    //! \brief access for sub entity data
    Data& operator () (const ElementType& element, const int subEntity ) 
    {
      return data_[ id_.subId( element, subEntity, codim_ ) ];
    }

    //! \brief access for sub entity data
    const Data& operator () (const ElementType& element, const int subEntity ) const 
    {
      return data_[ id_.subId( element, subEntity, codim_ ) ];
    }

    //! \brief iterator begin for iterating over data actually stored in container
    Iterator begin() 
    {
      return Iterator( data_.begin() );
    }

    //! \brief const iterator begin 
    ConstIterator begin() const
    {
      return ConstIterator( data_.begin() );
    }

    //! \brief iterator end  
    Iterator end() 
    {
      return Iterator( data_.end() );
    }

    //! \brief const iterator end  
    ConstIterator end() const
    {
      return ConstIterator( data_.end() );
    }

    //! \brief return size of allocated data 
    size_t size() const { return data_.size(); }

    //! \brief enlarge container, compress is not necessary but could be done
    void reserve()
    {
    }

    //! \brief adjust container to correct size and set all values to default
    void clear( )
    {
      data_.clear();
    }

    //! \brief adjust container to correct size including compress 
    void update( )
    { // this version could be implemented differently by only compressing
      update( Data() );
    }
  protected:  
    //! \brief adjust container to correct size including compress 
    void update( const Data& value )
    {
      // loop over all codimensions (needed to make codim_ static)
      ForLoop< AdaptCodim, 0, GridType :: dimension > :: apply( *this, value, codim_ );
    }

    template <int codim> 
    void adaptCodim( const Data& value )
    {
      assert( codim_ == codim );
      // create empty map and swap it with current map (no need to copy twice)
      Map oldData;
      std::swap( oldData, data_ );

      const iterator olddataend = oldData.end();
      typedef typename GridType :: template Codim< codim > :: LevelIterator LevelIterator ;
      typedef typename LevelIterator :: Entity  Entity; 
      for(int l = 0; l <= grid_.maxLevel(); ++ l) 
      {
        const LevelIterator endit = grid_.template lend< codim > ( l );   
        for( LevelIterator it = grid_.template lbegin< codim > ( l ); it != endit; ++ it )
        {
          const Entity& entity = * it ;
          const IdType id = id_.id( entity );
          Data& data = data_[ id ];
          iterator entry = oldData.find( id );
          if( entry != olddataend )
            data = (*entry).second;
        }
      }
    }
  };

  // PersistentContainer (default is to use PersistentContainerMap)
  // -------------------

  template < class Grid, class Data, class Allocator=std::allocator<Data> >  
  class PersistentContainer
  : public PersistentContainerMap< Grid, typename Grid::Traits::LocalIdSet, 
             std::map<const typename Grid::Traits::LocalIdSet::IdType, Data, 
                      std::less<const typename Grid::Traits::LocalIdSet::IdType>, Allocator> >
  {
    typedef Grid GridType;
    typedef typename Grid::Traits::LocalIdSet IdSet;
    typedef std::map<const typename Grid::Traits::LocalIdSet::IdType, Data,
                     std::less<const typename Grid::Traits::LocalIdSet::IdType>, Allocator> Map;
    typedef PersistentContainerMap< Grid, IdSet, Map > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim, grid.localIdSet() )
    {}
  };



  // PersistentContainer for AlbertaGrid
  // -----------------------------------

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< AlbertaGrid< dim, dimworld >, Data, Allocator >
  : public PersistentContainerVector< AlbertaGrid< dim, dimworld >,
                                      typename AlbertaGrid< dim, dimworld >::HierarchicIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef AlbertaGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1 )
    {}
  };



  // PersistentContainer for ALUGrid
  // -------------------------------

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUConformGrid< dim, dimworld >, Data, Allocator >
  : public PersistentContainerVector< ALUConformGrid< dim, dimworld >, 
                                      typename ALUConformGrid< dim, dimworld >::HierarchicIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef ALUConformGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1 )
    {}
  };

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUCubeGrid< dim, dimworld >, Data, Allocator >
  : public PersistentContainerVector< ALUCubeGrid< dim, dimworld >,
                                      typename ALUCubeGrid< dim, dimworld >::HierarchicIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef ALUCubeGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1 )
    {}
  };

  template< int dim, int dimworld, class Data, class Allocator >
  class PersistentContainer< ALUSimplexGrid< dim, dimworld >, Data, Allocator >
  : public PersistentContainerVector< ALUSimplexGrid< dim, dimworld >, 
                                      typename ALUSimplexGrid< dim, dimworld >::HierarchicIndexSet,
                                      std::vector<Data,Allocator> >
  {
    typedef ALUSimplexGrid< dim, dimworld > GridType;
    typedef PersistentContainerVector< GridType, typename GridType::HierarchicIndexSet, std::vector<Data,Allocator> > BaseType;

  public:
    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim, grid.hierarchicIndexSet(), 1.1 )
    {}
  };


  // PersistentContainer for GeometryGrid
  // ------------------------------------

#if 0 // not complete updated yet
  template< class HostGrid, class CoordFunction, class CoordAllocator, class Data, class Allocator >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction, CoordAllocator >, Data, Allocator >
  {
    typedef PersistentContainer< HostGrid, Data, Allocator > HostContainer;

  public:
    typedef GeometryGrid< HostGrid, CoordFunction, CoordAllocator > GridType;

    typedef typename HostContainer::ConstIterator ConstIterator;
    typedef typename HostContainer::Iterator Iterator;

    typedef typename GridType::template Codim< 0 >::Entity ElementType;

    //! Constructor filling the container with values using the default constructor 
    //! Depending on the implementation this could be achieved without allocating memory
    PersistentContainer ( const GridType &grid, const int codim )
    : hostContainer_( grid.hostGrid(), codim )
    {}
    //! Constructor filling the container with prescribed values;
    //! this will allocate memory for all entities of the given codimension
    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : hostContainer_( grid.hostGrid(), codim )
    {}

    template< class Entity >
    const Data &operator[] ( const Entity &entity ) const
    {
      return data( GridType::getRealImplementation( entity ) );
    }

    template< class Entity > 
    Data &operator[] ( const Entity &entity )
    {
      return data( GridType::getRealImplementation( entity ) );
    }

    const Data &operator() ( const ElementType &element, const int subEntity ) const
    {
      return hostContainer_( GridType::getRealImplementation( element ).hostEntity(), subEntity );
    }

    Data &operator() ( const ElementType &element, const int subEntity )
    {
      return hostContainer_( GridType::getRealImplementation( element ).hostEntity(), subEntity );
    }

    ConstIterator begin () const { return hostContainer_.begin(); }
    Iterator begin () { return hostContainer_.begin(); }

    ConstIterator end () const { return hostContainer_.end(); }
    Iterator end () { return hostContainer_.end(); }

    void resize ( const Data &value = Data() )
    {
      hostContainer_.resize( value );
    }

    void compress ( const Data &value = Data() )
    {
      hostContainer_.compress( value );
    }

  protected:
    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, Int2Type< false > ) const
    {
      return hostContainer_[ entity.hostEntity() ];
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, Int2Type< false > )
    {
      return hostContainer_[ entity.hostEntity() ];
    }

    template< class EntityImpl >
    const Data &data ( const EntityImpl &entity, Int2Type< true > ) const
    {
      return hostContainer_( entity.hostElement(). entity.subEntity() );
    }

    template< class EntityImpl >
    Data &data ( const EntityImpl &entity, Int2Type< true > )
    {
      return hostContainer_( entity.hostElement(). entity.subEntity() );
    }

  private:
    HostContainer &hostContainer_;
  };
#endif

} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
