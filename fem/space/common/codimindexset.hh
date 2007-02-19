#ifndef DUNE_CODIMINDEXSET_HH
#define DUNE_CODIMINDEXSET_HH

//- Dune includes 
#include <dune/common/misc.hh>
#include <dune/common/array.hh>

namespace Dune {
  
/*! 
 This class provides an Array which as an aditional resize and copy the old
 values functionality. 
*/
template <class T> 
class IndexArray : public Array<T> 
{
  int capacity_; 
public:
  IndexArray() : Array<T> () 
  {    
    this->n = 0; this->p = 0;
    capacity_ = this->n;
  }

  //! resize and set capacity 
  void resize(const int m) 
  {
    if( m <= capacity_ )
    {
      this->n = (m < 0) ? 0 : m;
      return;
    }
    Array<T> :: resize(m);
    capacity_ = this->n;
  }

  //! reallocate array with size m
  void realloc (const int m)
  {
    if(m <= this->n) 
    {
      // if m is smaller then zero, set size to zero 
      this->n = (m < 0) ? 0 : m;
      return;
    }
    
    int newSize = m*2;
    T * newMem = 0;
    
    try 
    {
      newMem = new T[newSize];
    }
    catch (std::bad_alloc) 
    {
      std::cerr << "Not enough memory!" << std::endl;
      throw;
    }

    if(this->n > 0)
      std::memcpy(newMem , this->p , this->n*sizeof(T));
    
    this->n = m;
    capacity_ = newSize;
    if(this->p) delete [] this->p;   
    this->p = newMem;
  }

  //! write Array to xdr stream
  bool processXdr(XDR *xdrs)
  {
    if(xdrs != 0)
    {
      int len = this->n;
      xdr_int( xdrs, &len );
      if(len != this->n) this->resize(len);
        xdr_vector(xdrs,(char *) this->p,this->n,sizeof(T),(xdrproc_t)xdr_double);
      return true;
    }
    else
      return false;
  }

};

//***********************************************************************
//
//  Index Set for one codimension
//  --CodimIndexSet 
//
//***********************************************************************
class CodimIndexSet
{
private:
  enum INDEXSTATE { NEW, USED, UNUSED };

  // the mapping of the global to leaf index 
  IndexArray<int> leafIndex_;

  // stack for holes 
  IndexArray<int> holes_; 
 
  // Array that only remeber the occuring holes (for compress of data)
  IndexArray<int> oldIdx_; 
  IndexArray<int> newIdx_; 

  // the state of each index 
  IndexArray<INDEXSTATE> state_;
 
  // next index to give away 
  int nextFreeIndex_;

  int actSize_;
  
  int myCodim_; 

  const double memFactor_;

public:
  //! Constructor
  CodimIndexSet (double memoryFactor = 1.1) 
    : nextFreeIndex_ (0)
    , actSize_(0)
    , myCodim_(-1) 
    , memFactor_(memoryFactor) 
  {
  }

  // set codim, because we can't use constructor 
  void setCodim (int codim) 
  {
    myCodim_ = codim;
  }

  // set codim, because we can't use constructor 
  int myCodim () const 
  {
    return myCodim_;
  }

  //! reallocate the vector for new size
  void resize ( const int newSize )
  {
    const int oldSize = leafIndex_.size();
    if(oldSize > newSize) return;

    const int newMemSize = ((int) memFactor_ * newSize);
    
    leafIndex_.realloc( newMemSize );
    state_.realloc( newMemSize );

    // reset new created parts of the vector 
    const int leafSize = leafIndex_.size();
    for(int i=oldSize; i<leafSize; ++i)
    {
      leafIndex_[i] = -1;
      state_[i] = UNUSED;
    }
  }

  //! set all entries to unused 
  void set2Unused() 
  {
    const int size = state_.size();
    for(int i=0; i<size; ++i) state_[i] = UNUSED;
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
    for(int i=0; i<sizeOfVecs; i++)
    {
      if(leafIndex_[i] >= 0)
      {
        // create vector with all holes 
        if(state_[i] == UNUSED)
        {
          holes_[actHole] = leafIndex_[i];
          actHole++;
        }

        // count the size of the leaf indices 
        newActSize++;
      }
    }

    assert( newActSize >= actHole );
    // the new size is the actual size minus the holes 
    actSize_ = newActSize - actHole;

    // resize hole storing vectors 
    oldIdx_.resize(actHole);
    newIdx_.resize(actHole);

    // close holes 
    //
    // NOTE: here the holes closing should be done in 
    // the opposite way. future work. 
     
    int holes = 0; // number of real holes 
    for(int i=0; i<leafIndex_.size(); i++)
    { 
      // a index that is used but larger then actual size 
      // has to move to a hole 
      if(state_[i] != UNUSED) 
      {
        // if used index lies behind size, then index has to move 
        // to one of the holes 
        if(leafIndex_[i] >= actSize_)
        {
          // serach next hole that is smaler than actual size 
          actHole--;
          // if actHole < 0 then error, because we have index larger then
          // actual size 
          assert(actHole >= 0);
          while ( holes_[actHole] >= actSize_ )
          {
            actHole--;
            if(actHole < 0) break;
          }

          assert(actHole >= 0);

#if HAVE_MPI 
          // only for none-ghost elements hole storage is applied
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
      else 
      {
        // all unsed indices are reset to -1 
        leafIndex_[i] = -1;
      }
    }

    // this call only sets the size of the vectors 
    oldIdx_.resize(holes);
    newIdx_.resize(holes);

    // the next index that can be given away is equal to size
    nextFreeIndex_ = actSize_;
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

  //! return leaf index for given hierarchic number  
  int index ( int num ) const
  {
    // assert if index was not set yet 
    return leafIndex_ [ num ];
  }
  
  //! return state of index for given hierarchic number  
  int state ( int num ) const
  {
    // assert if index was not set yet 
    assert( num >= 0 );
    assert( num < state_.size() );
    return state_ [ num ];
  }
  
  //! return state of index for given hierarchic number  
  bool exists ( int num ) const
  {
    return (state(num) != UNUSED);
  }
 
  //! return number of holes 
  int numberOfHoles () const
  {
    return oldIdx_.size();
  }

  //! return old index, for dof manager only 
  int oldIndex (int elNum ) const
  {
    return oldIdx_[elNum]; 
  }

  //! return new index, for dof manager only returns index 
  int newIndex (int elNum) const
  {
    return newIdx_[elNum]; 
  }

  // insert element and create index for element number  
  void insert (const int num )
  {
    assert(num < leafIndex_.size() );
    if(leafIndex_[num] < 0)
    {
      leafIndex_[num] = nextFreeIndex_;
      nextFreeIndex_++;
    }
    state_[num] = USED;
  }
  
  // insert element and create index for element number  
  void insertGhost (const int num )
  {
    assert(num < leafIndex_.size() );
    if(leafIndex_[num] < 0)
    {
      leafIndex_[num] = nextFreeIndex_;
      nextFreeIndex_++;
    }
    state_[num] = NEW;
  }
  
  // read/write from/to xdr stream 
  bool processXdr(XDR *xdrs)
  {
    xdr_int ( xdrs, &nextFreeIndex_ );
    xdr_int ( xdrs, &actSize_ );
    leafIndex_.processXdr(xdrs);
    state_.processXdr(xdrs);
    return true;
  }
  
  // insert element and create index for element number  
  void remove (int num )
  {
    assert(num < leafIndex_.size() );
    state_[num] = UNUSED;
  }
  
  void print (const char * msg, bool oldtoo = false ) const;

}; // end of CodimIndexSet  

} // end namespace Dune 
#endif
