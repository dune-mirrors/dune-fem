#ifndef DUNE_PERSISTENTCONTAINER_HH
#define DUNE_PERSISTENTCONTAINER_HH

#include <dune/fem/space/common/arrays.hh>

namespace Dune {

/** \brief class PersistentContainerVector */
template <class Grid, class Data>  
class PersistentContainerVector
{
protected:
  typedef typename Grid :: HierarchicIndexSet HierarchicIndexSetType;
  typedef Grid GridType;
  const GridType& grid_;
  const HierarchicIndexSetType& indexSet_;
  const int codim_;
  typedef std::vector< Data > StorageType;
  std::vector< Data > data_;
  
public:  
  typedef typename GridType :: template Codim< 0 > :: Entity ElementType; 
  typedef typename StorageType :: iterator Iterator ;
  typedef typename StorageType :: const_iterator ConstIterator ;

  //! constructor 
  PersistentContainerVector( const GridType& grid, const int codim , const Data& value = Data() ) 
    : grid_( grid )
    , indexSet_( grid_.hierarchicIndexSet() )
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
    : grid_( other.grid_ )
    , indexSet_( other.indexSet_ )
    , codim_( other.codim_ )
    , data_( other.data_ )  
  {}

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

  void set( const Data& value = Data() )
  {
    const size_t dataSize = data_.size();

    // set new value to default value 
    for(size_t i = 0; i<dataSize; ++i) 
      data_[ i ] = value;
  }

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

} // end namespace Dune   
#endif // end DUNE_PERSISTENTCONTAINER_HH
