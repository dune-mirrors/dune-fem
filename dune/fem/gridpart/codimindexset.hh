#ifndef DUNE_CODIMINDEXSET_HH
#define DUNE_CODIMINDEXSET_HH

#include <algorithm>

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

#include <dune/fem/io/streams/xdrstreams.hh>

#include <dune/fem/gridpart/persistentcontainer.hh>

namespace Dune {

//***********************************************************************
//
//  Index Set for one codimension
//  --CodimIndexSet 
//
//***********************************************************************
template <class GridImp>  
class CodimIndexSet
{
protected:  
  typedef GridImp GridType;
  typedef HierarchicIndexSetSelector< GridType > SelectorType;
  typedef typename SelectorType :: HierarchicIndexSet PersistentIndexSetType;

private:
  enum INDEXSTATE { UNUSED = 0,  // unused indices
                    USED   = 1,  // used indices
                    NEW    = 2 };//  new indices 

  // reference to persistent index container 
  const PersistentIndexSetType& indexContainer_;

  // array type for indices 
  typedef MutableArray<int> IndexArrayType;

  // Index pair that is stored, derive from pair to overload 
  // default constructor which sets the correct default data 
  struct IndexPair : public std::pair< int, INDEXSTATE > 
  {
    typedef std::pair< int, INDEXSTATE > BaseType;
    // default constructor 
    IndexPair() : BaseType( -1, UNUSED ) {}
  }; 

  typedef PersistentContainer< GridImp, IndexPair > IndexContainerType;

  // the mapping of the global to leaf index 
  IndexContainerType leafIndex_;

  // stack for holes 
  IndexArrayType holes_; 
 
  // Array that only remeber the occuring 
  // holes (for compress of data)
  IndexArrayType oldIdx_; 
  IndexArrayType newIdx_; 
 
  // next index to give away 
  int nextFreeIndex_;

  // last size of set before compress (needed in parallel runs) 
  int lastSize_;

  // codim for which index is provided 
  const int myCodim_; 

  // actual number of holes 
  int numberHoles_;

public:
  //! Constructor taking memory factor (default = 1.1)
  CodimIndexSet (const GridType& grid, 
                 const int codim, 
                 const double memoryFactor = 1.1) 
    : indexContainer_( SelectorType::hierarchicIndexSet( grid ) ) 
    , leafIndex_(grid, codim)
    , holes_(0)
    , oldIdx_(0)
    , newIdx_(0)
    , nextFreeIndex_ (0)
    , lastSize_ (0)
    , myCodim_( codim ) 
    , numberHoles_(0)
  {
    setMemoryFactor(memoryFactor);
  }

  //! set memory overestimation factor 
  void setMemoryFactor(const double memoryFactor)
  {
    holes_.setMemoryFactor(memoryFactor);
    oldIdx_.setMemoryFactor(memoryFactor);
    newIdx_.setMemoryFactor(memoryFactor);
  }

  //! returns vector with geometry tpyes this index set has indices for
  const std::vector <GeometryType> & geomTypes () const
  {
    return indexContainer_.geomTypes( myCodim_ );
  }

  //! reallocate the vectors
  void resize ()
  {
    leafIndex_.reserve();
  }

  //! prepare for setup (nothing to do here)
  void prepareCompress ()
  {
  }

public:  
  //! clear set 
  void clear() 
  {
    // set all values to -1 
    std::fill( leafIndex_.begin(), leafIndex_.end(), IndexPair() );
    // reset next free index 
    nextFreeIndex_ = 0;
  }

  //! set all entries to unused 
  void resetUsed() 
  {
    typedef typename IndexContainerType :: Iterator Iterator;
    const Iterator endit = leafIndex_.end();
    for( Iterator it = leafIndex_.begin(); it != endit; ++it )
    {
      (*it).second = UNUSED;
    }
  }

  //! set all entries to unused 
  void checkConsecutive() 
  {
    typedef typename IndexContainerType :: Iterator Iterator;
    const Iterator endit = leafIndex_.end();
    for( Iterator it = leafIndex_.begin(); it != endit; ++it )
    {
      const int idx = (*it).first; 
      assert( idx < nextFreeIndex_ );
    }
  }

  //! clear holes, i.e. set number of holes to zero 
  void clearHoles() 
  {
    // set number of holes to zero 
    numberHoles_ = 0;
    // remember actual size 
    lastSize_ = nextFreeIndex_;
  }

  //! make to index numbers consecutive 
  //! return true, if at least one hole was closed 
  bool compress ()
  {
    const int sizeOfVecs = leafIndex_.size();
    holes_.resize( sizeOfVecs );

    // true if a least one dof must be copied 
    bool haveToCopy = false;
    
    // mark holes 
    int actHole = 0;
    int newActSize = 0;
    typedef typename IndexContainerType :: Iterator Iterator;
    const Iterator endit = leafIndex_.end();
    for( Iterator it = leafIndex_.begin(); it != endit; ++it )
    {
      const IndexPair& leafIdx = *it;
      if( leafIdx.first >= 0 )
      {
        // create vector with all holes 
        if( leafIdx.second == UNUSED )
        {
          holes_[actHole] = leafIdx.first;
          ++actHole;
        }

        // count the size of the leaf indices 
        ++newActSize;
      }
    }

    assert( newActSize >= actHole );
    // the new size is the actual size minus the holes 
    int actSize = newActSize - actHole;

    // resize hole storing vectors 
    oldIdx_.resize(actHole);
    newIdx_.resize(actHole);

    // only compress if number of holes > 0    
    if(actHole > 0)
    {
      // close holes 
      //
      // NOTE: here the holes closing should be done in 
      // the opposite way. future work. 
      int holes = 0; // number of real holes 
      //size_t i = 0;
      const Iterator endit = leafIndex_.end();
      for( Iterator it = leafIndex_.begin(); it != endit; ++it )
      {
        IndexPair& leafIdx = *it;
        // a index that is used but larger then actual size 
        // has to move to a hole 
        if( leafIdx.second == UNUSED) 
        {
          // all unused indices are reset to -1 
          leafIdx.first = -1;
        }
        else 
        {
          // if used index lies behind size, then index has to move 
          // to one of the holes 
          if(leafIdx.first >= actSize)
          {
            // serach next hole that is smaler than actual size 
            actHole--;
            // if actHole < 0 then error, because we have index larger then
            // actual size 
            assert(actHole >= 0);
            while ( holes_[actHole] >= actSize )
            {
              actHole--;
              if(actHole < 0) break;
            }

            assert(actHole >= 0);

#if HAVE_MPI 
            // only for none-ghost elements hole storage is applied
            // this is because ghost indices might have in introduced 
            // after the resize was done. 
            if( leafIdx.second == USED ) 
#endif
            {
              // remember old and new index 
              oldIdx_[holes] = leafIdx.first; 
              newIdx_[holes] = holes_[actHole];
              ++holes;
            }
            
            leafIdx.first = holes_[actHole];

            // means that dof manager has to copy the mem
            leafIdx.second = NEW;
            haveToCopy = true;
          }
        }
      }

      // this call only sets the size of the vectors 
      oldIdx_.resize(holes);
      newIdx_.resize(holes);

    } // end if actHole > 0  
   
    // store number of actual holes 
    numberHoles_ = oldIdx_.size();

    // adjust size 
    leafIndex_.update();
    
    // the next index that can be given away is equal to size
    nextFreeIndex_ = actSize;

#ifndef NDEBUG
    checkConsecutive();
#endif

    return haveToCopy;
  }

  //! return how much extra memory is needed for restriction 
  int additionalSizeEstimate () const { return nextFreeIndex_; }

  //! return size of grid entities per level and codim 
  int size () const
  {
    return nextFreeIndex_;
  }
  
  //! return size of grid entities per level and codim 
  int realSize () const
  {
    return leafIndex_.size();
  }

  //! return leaf index for given entity   
  template <class EntityType>
  int index ( const EntityType& entity ) const
  {
    assert( myCodim_ == EntityType :: codimension );
    return leafIndex_[ entity ].first;
  }
  
  //! return leaf index for given entity   
  template <class EntityType>
  int subIndex ( const EntityType& entity,
                 const int subNumber ) const 
  {
    assert( 0 == EntityType :: codimension );
    return leafIndex_( entity, subNumber ).first;
  }
  
  //! return state of index for given hierarchic number  
  template <class EntityType> 
  bool exists ( const EntityType& entity ) const
  {
    assert( myCodim_ == EntityType :: codimension );
    return leafIndex_[ entity ].second != UNUSED;
  }
 
  //! return number of holes 
  int numberOfHoles () const
  {
    return numberHoles_;
  }

  //! return old index, for dof manager only 
  int oldIndex (int elNum ) const
  {
    assert( numberHoles_ == oldIdx_.size() );
    return oldIdx_[elNum]; 
  }

  //! return new index, for dof manager only returns index 
  int newIndex (int elNum) const
  {
    assert( numberHoles_ == newIdx_.size() );
    return newIdx_[elNum]; 
  }

  // insert element and create index for element number 
  template <class EntityType> 
  void insert (const EntityType& entity )
  {
    assert( myCodim_ == EntityType :: codimension );
    insertIdx( leafIndex_[ entity ] );
  }

  // insert element and create index for element number 
  template <class EntityType> 
  void insertSubEntity (const EntityType& entity, 
                        const int subNumber)  
  {
    assert( 0 == EntityType :: codimension );
    insertIdx( leafIndex_( entity, subNumber ) );
  }

  // insert element as ghost and create index for element number 
  template <class EntityType> 
  void insertGhost (const EntityType& entity )
  {
    assert( myCodim_ == EntityType :: codimension );
    // insert index 
    IndexPair& leafIdx = leafIndex_[ entity ];
    insertIdx( leafIdx );

    // if index is also larger than lastSize
    // mark as new to skip old-new index lists 
    if( leafIdx.first >= lastSize_ ) 
    {
      leafIdx.second = NEW;
    }
  }

  // insert element and create index for element number 
  template <class EntityType> 
  void markForRemoval( const EntityType& entity )
  {
    assert( myCodim_ == EntityType :: codimension );
    leafIndex_[ entity ].second = UNUSED;
  }

  // insert element as ghost and create index for element number 
  template <class EntityType> 
  bool validIndex (const EntityType& entity ) const
  {
    assert( myCodim_ == EntityType :: codimension );
    return leafIndex_[ entity ].first >= 0; 
  }

protected:
  // insert element and create index for element number  
  void insertIdx ( IndexPair& leafIdx )
  {
    if( leafIdx.first < 0 )
      leafIdx.first = nextFreeIndex_ ++ ;

    leafIdx.second = USED;
  }

public:  
  // write to stream 
  template <class StreamTraits> 
  bool write(OutStreamInterface< StreamTraits >& out) const
  {
    out << nextFreeIndex_;
    
    // backup leafIndex 
    //leafIndex_.write( out );

    return true;
  }
  
  // read from stream 
  template <class StreamTraits> 
  bool read(InStreamInterface< StreamTraits >& in)
  {
    in >> nextFreeIndex_;
    
    // restore leafIndex 
    //leafIndex_.read( in );

    return true;
  }
  
}; // end of CodimIndexSet  

} // end namespace Dune 
#endif
