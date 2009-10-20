#ifndef __OEM_TMPMEM_HH__
#define __OEM_TMPMEM_HH__

#include <cassert>

// temporay mem 
class OEMTmpMem
{
  double * mem_;
  int memSize_;
  int freeSize_;
public:
  OEMTmpMem () : mem_ (0) , memSize_ (0), freeSize_(0)
  {}

  ~OEMTmpMem ()
  {
    if(mem_) delete [] mem_;
  }

  void resize ( int newMemSize )
  {
    if(newMemSize > memSize_)
    {
      if(mem_) delete [] mem_;
      mem_ = new double [newMemSize];
      memSize_ = newMemSize;
      freeSize_ = memSize_;
    }
  }

  double * getMem ( int size )
  {
    freeSize_ -= size;
    assert(freeSize_ >= 0);
    return mem_ + freeSize_;
  }

  void reset ()
  {
    freeSize_ = memSize_;
  }
};

#endif

