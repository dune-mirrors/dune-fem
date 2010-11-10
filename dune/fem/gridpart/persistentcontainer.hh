#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <dune/fem/space/common/arrays.hh>

namespace Dune {

enum PersistentConainerComplexity { O_1 = 1, O_log_n = 2 };  

/** \brief class PersistentContainerVector */
template <class Grid, class Data>  
class PersistentContainerVector
{
protected:
  typedef typename Grid :: HierarchicIndexSet  IndexSetType;
  typedef Grid GridType;
  const IndexSetType& indexSet_;
  const int codim_;
  typedef std::vector< Data > StorageType;
  StorageType data_;
  
public:  
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 
  typedef typename StorageType :: iterator Iterator ;
  typedef typename StorageType :: const_iterator ConstIterator ;

  //! constructor 
  PersistentContainerVector( const GridType& grid, const int codim , const Data& value = Data() ) 
    : indexSet_( grid.hierarchicIndexSet() )
    , codim_( codim )
    , data_()
  {
    // resize to current size 
    adapt( value );
    // set memory overestimation factor 
    //data_.setMemoryFactor( 1.1 );
  }

  //! copy constructor 
  PersistentContainerVector( const PersistentContainerVector& other ) 
    : indexSet_( other.indexSet_ )
    , codim_( other.codim_ )
    , data_( other.data_ )  
  {}

  //! returns vector with geometry tpyes this index set has indices for
  const std::vector <GeometryType> & geomTypes () const
  {
    return indexSet_.geomTypes( codim_ );
  }
  
  PersistentConainerComplexity complexity() const { return O_1; }

  template <class Entity> 
  Data& operator [] (const Entity& entity ) 
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ indexSet_.index( entity ) ];
  }

  template <class Entity> 
  const Data& operator [] (const Entity& entity ) const
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ indexSet_.index( entity ) ];
  }

  Data& operator () (const ElementType& element, const int subEntity ) 
  {
    return data_[ indexSet_.subIndex( element, subEntity, codim_ ) ];
  }

  const Data& operator () (const ElementType& element, const int subEntity ) const 
  {
    return data_[ indexSet_.subIndex( element, subEntity, codim_ ) ];
  }

  size_t size() const { return data_.size(); }

  Iterator begin() 
  {
    return data_.begin();
  }

  ConstIterator begin() const
  {
    return data_.begin();
  }

  Iterator end() 
  {
    return data_.end();
  }

  ConstIterator end() const
  {
    return data_.end();
  }

  //! \brief adjsut container to new size 
  void resize(const Data& value = Data() )
  {
    if( indexSet_.size( codim_ ) > (int) data_.size() ) 
      adapt( value );
  }

  //! \brief adjsut container to new size 
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
  
public:  
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 
  typedef MyIterator< iterator > Iterator;
  typedef MyIterator< const_iterator > ConstIterator;

  //! constructor 
  PersistentContainerMap( const GridType& grid, const int codim , const Data& value = Data() ) 
    : grid_( grid )
    , idSet_( grid_.localIdSet() )
    , codim_( codim )
    , data_()
  {
    adapt( value );
  }

  //! copy constructor 
  PersistentContainerMap( const PersistentContainerMap& other ) 
    : grid_( other.grid_ )
    , idSet_( other.idSet_ )
    , codim_( other.codim_ )
    , data_( other.data_ )  
  {}

  PersistentConainerComplexity complexity() const { return O_log_n; }

  /*
  //! returns vector with geometry tpyes this index set has indices for
  const std::vector <GeometryType> & geomTypes () const
  {
    return indexContainer_.geomTypes( myCodim_ );
  }
  */
  
  template <class Entity> 
  Data& operator [] (const Entity& entity ) 
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ idSet_.id( entity ) ];
  }

  template <class Entity> 
  const Data& operator [] (const Entity& entity ) const
  { 
    assert( Entity :: codimension == codim_ );
    return data_[ idSet_.id( entity ) ];
  }

  Data& operator () (const ElementType& element, const int subEntity ) 
  {
    return data_[ idSet_.subId( element, subEntity, codim_ ) ];
  }

  const Data& operator () (const ElementType& element, const int subEntity ) const 
  {
    return data_[ idSet_.subId( element, subEntity, codim_ ) ];
  }

  size_t size() const { return data_.size(); }

  Iterator begin() 
  {
    return Iterator( data_.begin() );
  }

  ConstIterator begin() const
  {
    return ConstIterator( data_.begin() );
  }

  Iterator end() 
  {
    return Iterator( data_.end() );
  }

  ConstIterator end() const
  {
    return ConstIterator( data_.end() );
  }

  template <int codim> 
  void adaptCodim( const Data& value )
  {
    // make copy of data 
    StorageType oldData( data_ );

    // remove all data 
    data_.clear();

    typedef typename GridType :: template Codim< codim > :: LevelIterator LevelIterator ;
    typedef typename LevelIterator :: Entity  Entity; 
    const iterator olddataend = oldData.end();
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

  void resize( const Data& value = Data() )
  {
  }

  void adapt( const Data& value = Data() )
  {
    adaptCodim< 0 > ( value );
    // clear old data 
    // grid traversal 
  }
};

} // end namespace Dune   
#endif // end DUNE_PERSISTENTCONTAINER_HH
