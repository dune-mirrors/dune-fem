#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <map>
#include <vector>

#include <dune/common/misc.hh>
#include <dune/common/forloop.hh>
#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/common/arrays.hh>

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



enum PersistentConainerComplexity { O_1, O_log_n };  

/** \brief class PersistentContainerVector */
template <class Grid, class Data>  
class PersistentContainerVector
{
protected:
  typedef typename Grid :: HierarchicIndexSet  IndexSetType;
  typedef Grid GridType;
  const IndexSetType& indexSet_;
  const int codim_;
  //typedef std::vector< Data > StorageType;
  typedef MutableArray< Data > StorageType;
  StorageType data_;
  
public:  
  //! \brief entity of codimension 0 
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 

  //! \brief iterator type 
  typedef typename StorageType :: iterator Iterator ;
  //! \brief const iterator type 
  typedef typename StorageType :: const_iterator ConstIterator ;

  //! \brief constructor creating empty container 
  PersistentContainerVector( const GridType& grid, const int codim ) 
    : indexSet_( grid.hierarchicIndexSet() )
    , codim_( codim )
    , data_()
  {
    data_.setMemoryFactor( 1.1 );
  }

  //! \brief constructor also adapting to current grid size 
  PersistentContainerVector( const GridType& grid, const int codim, const Data& value ) 
    : indexSet_( grid.hierarchicIndexSet() )
    , codim_( codim )
    , data_()
  {
    data_.setMemoryFactor( 1.1 );
    adapt( value );
  }

  //! \brief copy constructor 
  PersistentContainerVector( const PersistentContainerVector& other ) 
    : indexSet_( other.indexSet_ )
    , codim_( other.codim_ )
    , data_( other.data_ )  
  {}

  //! \brief return complexity of random access 
  static PersistentConainerComplexity complexity () { return O_1; }

  //! \brief random access to entity data with correct codimension 
  template <class Entity> 
  Data& operator [] (const Entity& entity ) 
  { 
    assert( Entity :: codimension == codim_ );
    assert( indexSet_.index( entity ) < (typename IndexSetType::IndexType) data_.size() );
    return data_[ indexSet_.index( entity ) ];
  }

  //! \brief random access to entity data with correct codimension 
  template <class Entity> 
  const Data& operator [] (const Entity& entity ) const
  { 
    assert( Entity :: codimension == codim_ );
    assert( indexSet_.index( entity ) < (typename IndexSetType::IndexType) data_.size() );
    return data_[ indexSet_.index( entity ) ];
  }

  //! \brief access for sub entity data
  Data& operator () (const ElementType& element, const int subEntity ) 
  {
    assert( indexSet_.subIndex( element, subEntity, codim_ ) < (typename IndexSetType::IndexType) data_.size() );
    return data_[ indexSet_.subIndex( element, subEntity, codim_ ) ];
  }

  //! \brief access for sub entity data
  const Data& operator () (const ElementType& element, const int subEntity ) const 
  {
    assert( indexSet_.subIndex( element, subEntity, codim_ ) < (typename IndexSetType::IndexType) data_.size() );
    return data_[ indexSet_.subIndex( element, subEntity, codim_ ) ];
  }

  //! \brief return size of allocated data 
  size_t size() const { return data_.size(); }

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

  //! \brief enlarge container, compress is not necessary but could be done
  void enlarge( const Data &value = Data() )
  {
    if( (typename IndexSetType::IndexType) indexSet_.size( codim_ ) > (typename IndexSetType::IndexType) data_.size() ) 
      adapt( value );
  }

  //! \brief adjust container to correct size including compress 
  void adapt(const Data& value = Data() )
  {
    const size_t oldSize = data_.size();
    const size_t dataSize = indexSet_.size( codim_ );
    data_.resize( dataSize );

    // set new value to default value 
    for(size_t i = oldSize; i<dataSize; ++i) 
      data_[ i ] = value;
  }
};

/** \brief class PersistentContainerVector */
template <class Grid, class Data>  
class PersistentContainerMap
{
  typedef PersistentContainerMap< Grid, Data > ThisType;

protected:
  typedef typename Grid :: Traits :: LocalIdSet IdSetType;
  typedef typename IdSetType :: IdType  IdType;
  typedef Grid GridType;

  const GridType& grid_;
  const IdSetType& idSet_;
  const int codim_;
  typedef std::map< const IdType, Data > StorageType;
  mutable StorageType data_;

  typedef typename StorageType :: iterator iterator ;
  typedef typename StorageType :: const_iterator const_iterator ;

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

  //! \brief constructor creating empty container 
  PersistentContainerMap( const GridType& grid, const int codim ) 
    : grid_( grid )
    , idSet_( grid_.localIdSet() )
    , codim_( codim )
    , data_()
  {
  }

  //! \brief constructor also adapting to current grid size 
  PersistentContainerMap( const GridType& grid, const int codim , const Data& value ) 
    : grid_( grid )
    , idSet_( grid_.localIdSet() )
    , codim_( codim )
    , data_()
  {
    adapt( value );
  }

  //! \brief copy constructor 
  PersistentContainerMap( const PersistentContainerMap& other ) 
    : grid_( other.grid_ )
    , idSet_( other.idSet_ )
    , codim_( other.codim_ )
    , data_( other.data_ )  
  {}

  //! \brief return complexity of random access 
  static PersistentConainerComplexity complexity () { return O_log_n; }

  //! \brief random access entity with correct codimension 
  template <class Entity> 
  Data& operator [] (const Entity& entity ) 
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ idSet_.id( entity ) ];
  }

  //! \brief random access entity with correct codimension 
  template <class Entity> 
  const Data& operator [] (const Entity& entity ) const
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ idSet_.id( entity ) ];
  }

  //! \brief access for sub entity data
  Data& operator () (const ElementType& element, const int subEntity ) 
  {
    return data_[ idSet_.subId( element, subEntity, codim_ ) ];
  }

  //! \brief access for sub entity data
  const Data& operator () (const ElementType& element, const int subEntity ) const 
  {
    return data_[ idSet_.subId( element, subEntity, codim_ ) ];
  }

  //! \brief return size of allocated data 
  size_t size() const { return data_.size(); }

  //! \brief iterator begin 
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

  //! \brief enlarge container, compress is not necessary but could be done
  void enlarge( const Data& value = Data() )
  {
  }

  //! \brief adjust container to correct size including compress 
  void adapt( const Data& value = Data() )
  {
    // loop over all codimensions 
    ForLoop< AdaptCodim, 0, GridType :: dimension > :: apply( *this, value, codim_ );
  }

protected:  
  template <int codim> 
  void adaptCodim( const Data& value )
  {
    assert( codim_ == codim );
    // create empty map and swap it with current map (no need to copy twice)
    StorageType oldData;
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
        const IdType id = idSet_.id( entity );
        Data& data = data_[ id ];
        iterator entry = oldData.find( id );
        if( entry == olddataend )
          data = value ;
        else 
          data = (*entry).second;
      }
    }
  }
};



  // PersistentContainer
  // -------------------

  template< class Grid, class Data >
  class PersistentContainer
  : public PersistentContainerMap< Grid, Data >
  {
    typedef PersistentContainerMap< Grid, Data > BaseType;

  public:
    PersistentContainer ( const Grid &grid, const int codim )
    : BaseType( grid, codim )
    {}

    PersistentContainer ( const Grid &grid, const int codim, const Data &value )
    : BaseType( grid, codim, value )
    {}
  };



  // PersistentContainer for AlbertaGrid
  // -----------------------------------

  template< int dim, int dimworld, class Data >
  class PersistentContainer< AlbertaGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< AlbertaGrid< dim, dimworld >, Data >
  {
    typedef PersistentContainerVector< AlbertaGrid< dim, dimworld >, Data > BaseType;

  public:
    typedef AlbertaGrid< dim, dimworld > GridType;

    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim )
    {}

    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : BaseType( grid, codim, value )
    {}
  };



  // PersistentContainer for ALUGrid
  // -------------------------------

  template< int dim, int dimworld, class Data >
  class PersistentContainer< ALUConformGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< ALUConformGrid< dim, dimworld >, Data >
  {
    typedef PersistentContainerVector< ALUConformGrid< dim, dimworld >, Data > BaseType;

  public:
    typedef ALUConformGrid< dim, dimworld > GridType;

    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim )
    {}

    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : BaseType( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, class Data >
  class PersistentContainer< ALUCubeGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< ALUCubeGrid< dim, dimworld >, Data >
  {
    typedef PersistentContainerVector< ALUCubeGrid< dim, dimworld >, Data > BaseType;

  public:
    typedef ALUCubeGrid< dim, dimworld > GridType;

    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim )
    {}

    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : BaseType( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, class Data >
  class PersistentContainer< ALUSimplexGrid< dim, dimworld >, Data >
  : public PersistentContainerVector< ALUSimplexGrid< dim, dimworld >, Data >
  {
    typedef PersistentContainerVector< ALUSimplexGrid< dim, dimworld >, Data > BaseType;

  public:
    typedef ALUSimplexGrid< dim, dimworld > GridType;

    PersistentContainer ( const GridType &grid, const int codim )
    : BaseType( grid, codim )
    {}

    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : BaseType( grid, codim, value )
    {}
  };



  // PersistentContainer for GeometryGrid
  // ------------------------------------

  template< class HostGrid, class CoordFunction, class Allocator, class Data >
  class PersistentContainer< GeometryGrid< HostGrid, CoordFunction, Allocator >, Data >
  {
    typedef PersistentContainer< HostGrid, Data > HostContainer;

  public:
    typedef GeometryGrid< HostGrid, CoordFunction, Allocator > GridType;

    typedef typename HostContainer::ConstIterator ConstIterator;
    typedef typename HostContainer::Iterator Iterator;

    typedef typename GridType::template Codim< 0 >::Entity ElementType;

    PersistentContainer ( const GridType &grid, const int codim )
    : hostContainer_( grid.hostGrid(), codim )
    {}

    PersistentContainer ( const GridType &grid, const int codim, const Data &value )
    : hostContainer_( grid.hostGrid(), codim, value )
    {}

    static PersistentConainerComplexity complexity ()
    {
      return HostContainer::complexity();
    }

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

    size_t size () const { return hostContainer_.size(); }

    ConstIterator begin () const { return hostContainer_.begin(); }
    Iterator begin () { return hostContainer_.begin(); }

    ConstIterator end () const { return hostContainer_.end(); }
    Iterator end () { return hostContainer_.end(); }

    void enlarge ( const Data &value = Data() )
    {
      hostContainer_.enlarge( value );
    }

    void adapt ( const Data &value = Data() )
    {
      hostContainer_.adapt( value );
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

} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
