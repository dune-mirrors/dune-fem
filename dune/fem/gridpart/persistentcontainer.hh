#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <map>
#include <vector>

#include <dune/common/forloop.hh>
#include <dune/grid/common/capabilities.hh>

#include <dune/fem/space/common/arrays.hh>

namespace Dune {

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
    if( indexSet_.size( codim_ ) > (typename IndexSetType::IndexType) data_.size() ) 
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

} // end namespace Dune

#endif // end DUNE_PERSISTENTCONTAINER_HH
