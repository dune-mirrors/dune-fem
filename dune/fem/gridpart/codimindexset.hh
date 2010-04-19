#ifndef DUNE_CODIMINDEXSET_HH
#define DUNE_CODIMINDEXSET_HH

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/gridpart/defaultindexsets.hh>

#if DUNE_FEM_COMPATIBILITY
#include <dune/fem/io/file/xdrio.hh>
#endif

#include <dune/fem/io/streams/xdrstreams.hh>

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
  //enum INDEXSTATE { NEW = 2 , USED = 1 , UNUSED = -1 };
  static const char NEW    = 2;  // new indices 
  static const char USED   = 1;  // used indices 
  static const char UNUSED = 0; // unused indices 

  // reference to persistent index container 
  const PersistentIndexSetType& indexContainer_;

  // array type for indices 
  typedef MutableArray<int> IndexArrayType;

  // array type of state of indices 
  typedef MutableArray<char> StateArrayType;

  // the mapping of the global to leaf index 
  IndexArrayType leafIndex_;

  // stack for holes 
  IndexArrayType holes_; 
 
  // Array that only remeber the occuring 
  // holes (for compress of data)
  IndexArrayType oldIdx_; 
  IndexArrayType newIdx_; 

  // the state of each index 
  StateArrayType state_;
 
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
    , leafIndex_(0)
    , holes_(0)
    , oldIdx_(0)
    , newIdx_(0)
    , state_(0)
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
    leafIndex_.setMemoryFactor(memoryFactor);
    state_.setMemoryFactor(memoryFactor);
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
    resize( indexContainer_.size( myCodim_ ) );
  }

protected: 
  //! reallocate the vector for new size
  void resize ( const int newSize )
  {
    const int oldSize = leafIndex_.size();
    if(oldSize > newSize) return;

    leafIndex_.resize( newSize );
    state_.resize( newSize );

    // reset new created parts of the vector 
    const int leafSize = leafIndex_.size();
    for(int i=oldSize; i<leafSize; ++i)
    {
      leafIndex_[i] = -1;
      state_[i] = UNUSED;
    }
  }

public:  
  //! clear set 
  void clear() 
  {
    // remove all information 
    {
      const int size = state_.size();
      for(int i=0; i<size; ++i) state_[i] = UNUSED;
    }

    {
      const int size = leafIndex_.size();
      for(int i=0; i<size; ++i) leafIndex_[i] = -1;
    }
    // reset next free index 
    nextFreeIndex_ = 0;
  }

  //! set all entries to unused 
  void resetUsed() 
  {
    const int size = state_.size();
    for(int i=0; i<size; ++i) state_[i] = UNUSED;
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
    const int sizeOfVecs = state_.size();
    holes_.resize( sizeOfVecs );

    // true if a least one dof must be copied 
    bool haveToCopy = false;
    
    // mark holes 
    int actHole = 0;
    int newActSize = 0;
    for(int i=0; i<sizeOfVecs; ++i)
    {
      if(leafIndex_[i] >= 0)
      {
        // create vector with all holes 
        if(state_[i] == UNUSED)
        {
          holes_[actHole] = leafIndex_[i];
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
      for(int i=0; i<leafIndex_.size(); ++i)
      { 
        // a index that is used but larger then actual size 
        // has to move to a hole 
        if(state_[i] == UNUSED) 
        {
          // all unused indices are reset to -1 
          leafIndex_[i] = -1;
        }
        else 
        {
          // if used index lies behind size, then index has to move 
          // to one of the holes 
          if(leafIndex_[i] >= actSize)
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
            if( state_[i] == USED ) 
#endif
            {
              // remember old and new index 
              oldIdx_[holes] = leafIndex_[i]; 
              newIdx_[holes] = holes_[actHole];
              ++holes;
            }
            
            leafIndex_[i] = holes_[actHole];

            // means that dof manager has to copy the mem
            state_[i] = NEW;
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

    // the next index that can be given away is equal to size
    nextFreeIndex_ = actSize;
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
    return indexIdx( indexContainer_.index( entity ) );
  }
  
  //! return leaf index for given entity   
  template <class EntityType>
  int subIndex ( const EntityType& entity,
                 const int subNumber ) const 
  {
    assert( 0 == EntityType :: codimension );
    return indexIdx( indexContainer_.subIndex( entity, subNumber, myCodim_ ) );
  }
  
  //! return state of index for given hierarchic number  
  template <class EntityType> 
  bool exists ( const EntityType& entity ) const
  {
    assert( myCodim_ == EntityType :: codimension );
    return existsIdx( indexContainer_.index( entity ) );
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
    return insertIdx( indexContainer_.index( entity ) );
  }

  // insert element and create index for element number 
  template <class EntityType> 
  void insertSubEntity (const EntityType& entity, 
                        const int subNumber)  
  {
    assert( 0 == EntityType :: codimension );
    insertIdx( indexContainer_.subIndex( entity, subNumber, myCodim_ ) );
  }

  // insert element as ghost and create index for element number 
  template <class EntityType> 
  void insertGhost (const EntityType& entity )
  {
    assert( myCodim_ == EntityType :: codimension );
    insertGhostIdx( indexContainer_.index( entity ) );
  }

  // insert element and create index for element number 
  template <class EntityType> 
  void remove( const EntityType& entity )
  {
    assert( myCodim_ == EntityType :: codimension );
    removeIdx( indexContainer_.index( entity ) );
  }

  // insert element as ghost and create index for element number 
  template <class EntityType> 
  bool validIndex (const EntityType& entity ) const
  {
    assert( myCodim_ == EntityType :: codimension );
    return validIndexIdx( indexContainer_.index( entity ) );
  }

protected:
  //! return true if index is valid 
  bool validIndexIdx ( const int num ) const
  {
    return (leafIndex_[num] >= 0);
  }
 
  //! return leaf index for given hierarchic number  
  int indexIdx ( const int num ) const
  {
    // assert if index was not set yet 
    return leafIndex_ [ num ];
  }
  
  //! return state of index for given hierarchic number  
  bool existsIdx ( const int num ) const
  {
    // assert if index was not set yet 
    assert( num >= 0 );
    assert( num < state_.size() );
    return state_ [ num ] != UNUSED;
  }
 
  // insert element and create index for element number  
  void insertIdx (const int num )
  {
    assert(num < leafIndex_.size() );
    if(leafIndex_[num] < 0)
    {
      leafIndex_[num] = nextFreeIndex_;
      ++nextFreeIndex_;
    }
    state_[num] = USED;
  }
  
  // insert element as ghost and create index for element number  
  void insertGhostIdx (const int num )
  {
    // insert index 
    insertIdx( num );

    // if index is also larger than lastSize
    // mark as new to skip old-new index lists 
    if( leafIndex_[num] >= lastSize_ ) 
    {
      state_[num] = NEW;
    }
  }
  
  // insert element and create index for element number  
  void removeIdx (const int num)
  {
    assert(num < leafIndex_.size() );
    state_[num] = UNUSED;
  }
public:  
#if DUNE_FEM_COMPATIBILITY
  // read/write from/to xdr stream 
  bool processXdr(XDRStream& xdr)
  {
    // restore size of index set 
    int ret = xdr.inout( nextFreeIndex_ );
    bool ok = (ret == 1) ? true : false;
    
    // should always have the same length 
    assert( leafIndex_.size() == state_.size() );

    // backup/restore leafIndex 
    ok |= leafIndex_.processXdr(xdr);
    // backup/restore state 
    ok |= state_.processXdr(xdr);

    // should always have the same length 
    assert( leafIndex_.size() == state_.size() );

    return ok;
  }
#endif
  
  // write to stream 
  template <class StreamTraits> 
  bool write(OutStreamInterface< StreamTraits >& out) const
  {
    out << nextFreeIndex_;
    
    // should always have the same length 
    assert( leafIndex_.size() == state_.size() );

    // backup leafIndex 
    leafIndex_.write( out );

    // backup state
    state_.write( out );

    return true;
  }
  
  // read from stream 
  template <class StreamTraits> 
  bool read(InStreamInterface< StreamTraits >& in)
  {
    in >> nextFreeIndex_;
    
    // restore leafIndex 
    leafIndex_.read( in );

    // restore state
    state_.read( in );

    // should always have the same length 
    assert( leafIndex_.size() == state_.size() );

    return true;
  }
  
}; // end of CodimIndexSet  

} // end namespace Dune 
#endif
