#ifndef DUNE_IDBASEDCODIMINDEXSET_HH
#define DUNE_IDBASEDCODIMINDEXSET_HH

//- system includes 
#include <map>
#include <stack>
#include <vector>

#include <dune/common/geometrytype.hh>
#include <dune/fem/space/common/arrays.hh>

namespace Dune
{

  //! Send vector to output stream
  template<class MapType>
  inline void printMap (std::ostream& s, const MapType& m)
  {
    typedef typename MapType :: const_iterator iterator;
    iterator end = m.end();
    for(iterator it=m.begin(); it != end; ++it) 
    {
      int idx = (*it).second;
      s << idx << " ";  
    }
    s << std::endl;
  }



  //! \brief LeafIndexSet based on local ids using std::map to build mapping
  //!from id to index 
  template< class GridImp >
  class IdBasedCodimIndexSet
  {
    typedef IdBasedCodimIndexSet< GridImp > ThisType;
  protected:  
    typedef GridImp  GridType;
    // use LocalIdSet to get persistent ids 
    typedef typename GridType :: Traits :: LocalIdSet   IdSetType;
    typedef typename GridType :: Traits :: LeafIndexSet IndexSetType;

    typedef typename IdSetType::IdType IdType;
    
    typedef std::map< IdType, int  > IndexStorageType;
    typedef std::map< int , IdType > IdStorageType;

    typedef MutableArray< int > IndexVectorType;
    
    // the mapping of the global to leaf index 
    mutable IndexStorageType* leafIndex_;

    // old indices 
    mutable IndexStorageType* oldLeafIndex_;

    // list of old and new index for all holes 
    mutable IndexVectorType oldIndexVec_;
    mutable IndexVectorType newIndexVec_;

    const GridType     & grid_;
    const IdSetType    & idSet_;
    const IndexSetType & indexSet_;

    // my codimension 
    const int codim_; 
   
    // next index to give away 
    int nextFreeIndex_;

    // number of holes 
    int numberOfHoles_; 

    // no copying 
    IdBasedCodimIndexSet (const IdBasedCodimIndexSet& );

  public:
    //! Constructor
    IdBasedCodimIndexSet ( const GridType& grid,
                           const int codim, 
                           const double memFactor = 1.0 )
    : leafIndex_( new IndexStorageType () ),
      oldLeafIndex_( new IndexStorageType () ),
      grid_( grid ),
      idSet_( grid.localIdSet() ),
      indexSet_( grid.leafIndexSet() ),
      codim_( codim ),
      nextFreeIndex_( 0 )
    {
      // set memory over estimation 
      oldIndexVec_.setMemoryFactor( 1.1 );
      newIndexVec_.setMemoryFactor( 1.1 );
    }

    //! destructor 
    ~IdBasedCodimIndexSet() 
    {
      delete oldLeafIndex_ ;
      delete leafIndex_ ;
    }

    //! set all entries to unused 
    void resetUsed() {}

    //! returns vector with geometry tpyes this index set has indices for
    const std::vector <GeometryType> & geomTypes () const
    {
      return indexSet_.geomTypes( codim_ );
    }

    //! nothing to do here
    void resize ()
    {
    }

    void prepareCompress()
    {
      // swap pointers of index storages 
      // this way we keep the old pattern and can 
      // simply fill with new indices 
      {
        IndexStorageType* swap  = oldLeafIndex_;
        oldLeafIndex_ = leafIndex_;
        leafIndex_ = swap;
      }

      leafIndex().clear();
    }

    //! make to index numbers consecutive 
    //! return true, if at least one hole was closed 
    bool compress ()
    {
      const int currentSize = indexSet_.size( codim_ ); 
      nextFreeIndex_ = currentSize;

      // create holes index vector 
      IndexVectorType holesIdx(leafIndex().size());   

      int noHole = 0;
      {
        // check for old indices that are smaller than current size 
        typedef typename IndexStorageType :: const_iterator iterator;
        iterator end = oldLeafIndex().end();
        for(iterator it = oldLeafIndex().begin(); it != end; ++it) 
        {
          const int idx = (*it).second;
          assert( idx >= 0 );
          
          // check if index is larger than currentSize it is a hole 
          if( (idx < currentSize) && (idx >= 0) )
          {
            holesIdx[ noHole ] = idx;
            ++noHole;
          }
        }
      }

      {
        // resize hole lists 
        oldIndexVec_.resize( leafIndex().size() );
        newIndexVec_.resize( leafIndex().size() );
        
        int hole = 0;
        typedef typename IndexStorageType :: iterator iterator;
        const iterator end = leafIndex().end();
        for(iterator it = leafIndex().begin(); it != end; ++it) 
        {
          const int idx = (*it).second;
          if( idx >= nextFreeIndex_ ) 
          {
            assert( hole < noHole );
            assert( hole < oldIndexVec_.size() );
            assert( hole < newIndexVec_.size() );

            // put to old idx list 
            oldIndexVec_[ hole ] = idx;
            
            const int newIdx = holesIdx[ hole ];
            
            // store new index 
            newIndexVec_[ hole ] = newIdx;
            
            // set new index to index map 
            (*it).second = newIdx;

            // next hole 
            ++hole;
          }

          assert( (*it).second < size() ); 
        }

        // set current number of holes
        oldIndexVec_.resize( hole );
        newIndexVec_.resize( hole );

        // store number of holes 
        numberOfHoles_ = oldIndexVec_.size();
      }

      // clear old values 
      oldLeafIndex().clear();
      
      // check that index set is consecutive 
      // only done in debug mode 
      checkConsecutive();
      return true; 
    }

    void checkConsecutive()
    {
  #ifndef NDEBUG
      typedef typename IndexStorageType :: const_iterator iterator;
      iterator end    = leafIndex().end();
      for(iterator it = leafIndex().begin(); it != end; ++it) 
      {
        const int idx = (*it).second;
        assert( idx < nextFreeIndex_ );
      }
  #endif
    }

    //! return how much extra memory is needed for restriction 
    int additionalSizeEstimate () const { return size(); }

    //! return size of grid entities per level and codim 
    int size () const
    {
      return nextFreeIndex_;
    }
    
    //! return size of grid entities per level and codim 
    int realSize () const
    {
      return size();
    }

    //! return leaf index for given hierarchic number  
    int indexId ( const IdType& id ) const
    {
      // assert if index was not set yet 
      assert( leafIndex() [ id ] < size() );
      return leafIndex() [ id ];
    }
   
    //! return leaf index for given hierarchic number  
    template <class EntityType> 
    int index ( const EntityType & en ) const
    {
      // assert if index was not set yet 
      assert( EntityType :: codimension == codim_ );
      assert( exists( en ) );
      return indexId( idSet_.id(en) );
    }
   
    //! return leaf index for given hierarchic number  
    template <class EntityType> 
    int subIndex ( const EntityType & en, const int subNumber) const
    {
      assert( EntityType :: codimension == 0 );
      return indexId( idSet_.subId( en, subNumber, codim_ ) );
    }
   
    //! return state of index for given hierarchic number  
    template <class EntityType> 
    bool exists ( const EntityType & en) const
    {
      assert( codim_ == EntityType :: codimension );
      return (leafIndex().find( idSet_.id( en ) ) != leafIndex().end()); 
    }

    // return true if index of entity is valid 
    template <class EntityType>
    bool validIndex (const EntityType& entity ) const
    {
      return exists( entity );
    }

    void clear() {}
   
    //! remove the list of holes (to avoid double compression)
    void clearHoles() 
    {
      numberOfHoles_ = 0;
    }

    //! return number of existing holes 
    int numberOfHoles () const
    {
      assert( oldIndexVec_.size() == newIndexVec_.size() );
      return numberOfHoles_; 
    }

    //! return old index, for dof manager only 
    int oldIndex ( const int idx ) const
    {
      return oldIndexVec_[ idx ];
    }

    //! return new index, for dof manager only returns index 
    int newIndex ( const int idx ) const
    {
      return newIndexVec_[ idx ];
    }
    
    // check index   
    template <class EntityType> 
    void insert(const EntityType & en) 
    {
      insertId( idSet_.id(en) );
    }

    // insert element and create index for element number  
    void insertId ( const IdType& id ) 
    {
      if( leafIndex().find(id) == leafIndex().end() )
      {
        // check if entity is really new 
        if( oldLeafIndex().find(id) == oldLeafIndex().end() )
        {
          const int idx = nextFreeIndex_;
          ++nextFreeIndex_;
          leafIndex()[ id  ] = idx; 
        }
        else 
        {
          const int idx = oldLeafIndex()[ id ];
          // remove from list of unused indices 
          oldLeafIndex().erase( id );
          leafIndex()[ id  ] = idx; 
        }
      }
    }
    
    // insert element and create index for element number  
    template <class EntityType> 
    void insertGhost (const EntityType & en ) 
    {
      insert( en );
    }

    // insert element and create index for element number 
    template <class EntityType>
    void insertSubEntity (const EntityType& entity,
                          const int subNumber)
    {
      assert( 0 == EntityType :: codimension );
      insertId( idSet_.subId( entity, subNumber, codim_ ) );
    }
    
    // removal of indices is actually done in compression 
    template <class EntityType> 
    void markForRemoval( const EntityType & en ) 
    {
    }
    
    // write to stream 
    template <class StreamTraits>
    bool write(OutStreamInterface< StreamTraits >& out) const
    {
      return true;
    }

    // read from stream 
    template <class StreamTraits>
    bool read(InStreamInterface< StreamTraits >& in)
    {
      return true;
    }
      
  protected:
    // return reference to leaf index 
    IndexStorageType& leafIndex() const 
    {
      assert( leafIndex_ );
      return *leafIndex_;
    }

    // return reference to old leaf index 
    IndexStorageType& oldLeafIndex() const 
    {
      assert( oldLeafIndex_ );
      return *oldLeafIndex_;
    }

  }; // end of class 

} // end namespace Dune 

#endif
