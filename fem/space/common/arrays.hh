#ifndef DUNE_ARRAYS_HH
#define DUNE_ARRAYS_HH

//- System includes 
#include <cassert>
#include <string>

//- Dune includes 
#include <dune/common/genericiterator.hh>
#include <dune/common/interfaces.hh>

#if HAVE_BLAS 
// include BLAS for daxpy operation 
#include <dune/fem/solver/oemsolver/cblas.h>
#endif

namespace Dune {

// forward declarations 
template <class T>
class DefaultDofAllocator;
  
template <class T, class AllocatorType = DefaultDofAllocator<T> >
class MutableArray;

//! oriented to the STL Allocator funtionality 
template <class T>
class DefaultDofAllocator {
public:
  //! allocate array of nmemb objects of type T
  static T* malloc (size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = new T[nmemb];
    assert( p );
    return p;
  }

  //! release memory previously allocated with malloc member
  static void free (T* p)
  {
    delete [] p;
  }
  
  //! allocate array of nmemb objects of type T
  static T* realloc (T* oldMem, size_t oldSize , size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = DefaultDofAllocator<T> ::malloc(nmemb);
    std::memcpy(p, oldMem, oldSize * sizeof(T) );
    DefaultDofAllocator<T> :: free (oldMem);
    return p;
  }
};

//! allocator for simple structures like int, double and float
//! using the C malloc,free, and realloc 
struct SimpleDofAllocator 
{
  //! allocate array of nmemb objects of type T
  template <typename T>
  static T* malloc (size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = (T *) std::malloc(nmemb * sizeof(T));
    assert(p);
    return p;
  }

  //! release memory previously allocated with malloc member
  template <typename T>
  static void free (T* p)
  {
    assert(p);
    std::free(p);
  }
  
  //! allocate array of nmemb objects of type T
  template <typename T>
  static T* realloc (T* oldMem, size_t oldSize , size_t nmemb)
  {
    assert(nmemb > 0);
    assert(oldMem);
    T * p = (T *) std::realloc(oldMem , nmemb*sizeof(T));
    assert(p);
    return p;
  }
};

template <>
class DefaultDofAllocator<double> 
{
  typedef double T;
public:
  //! allocate array of nmemb objects of type T
  static T* malloc (size_t nmemb)
  {
    return SimpleDofAllocator::malloc<T> (nmemb);
  }

  //! release memory previously allocated with malloc member
  static void free (T* p)
  {
    SimpleDofAllocator::free<T> (p);
    return ;
  }
  
  //! allocate array of nmemb objects of type T
  static T* realloc (T* oldMem, size_t oldSize , size_t nmemb)
  {
    return SimpleDofAllocator::realloc<T> (oldMem,oldSize,nmemb);
  }
};

template <>
class DefaultDofAllocator<int> 
{
  typedef int T;
public:
  //! allocate array of nmemb objects of type T
  static T* malloc (size_t nmemb)
  {
    return SimpleDofAllocator::malloc<T> (nmemb);
  }

  //! release memory previously allocated with malloc member
  static void free (T* p)
  {
    SimpleDofAllocator::free<T> (p);
    return ;
  }
  
  //! allocate array of nmemb objects of type T
  static T* realloc (T* oldMem, size_t oldSize , size_t nmemb)
  {
    return SimpleDofAllocator::realloc<T> (oldMem,oldSize,nmemb);
  }
};

//! specialisation for MutableArray<V> that uses std::copy for copying
//! memory 
template <class V>
class DefaultDofAllocator<MutableArray<V> > 
{
  typedef MutableArray<V> T;
public:
  //! allocate array of nmemb objects of type T
  static T* malloc (size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = new T[nmemb];
    assert( p );
    return p;
  }

  //! release memory previously allocated with malloc member
  static void free (T* p)
  {
    delete [] p;
  }
  
  //! allocate array of nmemb objects of type T
  static T* realloc (T* oldMem, size_t oldSize , size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = DefaultDofAllocator<T> ::malloc(nmemb);
    std::copy( oldMem->begin(), oldMem->end(), p->begin());
    DefaultDofAllocator<T> :: free (oldMem);
    return p;
  }
};

/** \brief Static Array Wrapper for simple C Vectors like double* and
  int*. This also works as base class for the MutableArray which is used
  to store the degrees of freedom. 
*/
template <class T> 
class StaticArray
{
protected:
  typedef StaticArray<T> ThisType;

  // size of array 
  int size_;

  // pointer to mem
  T * vec_;

  StaticArray(const StaticArray&);
public:
  //! definition conforming to STL  
  typedef T value_type;
  
  //! DofIterator
  typedef GenericIterator<ThisType, T> DofIteratorType;
  
  //! Const DofIterator
  typedef GenericIterator<const ThisType, const T> ConstDofIteratorType;

  //! create array of length size
  //! if size is <= 0 then vec of lenght 1 is created (parallel runs)
  StaticArray(const int size, T* vec) 
    : size_(size)
    , vec_(vec) 
  {
    assert( size_ >= 0 );
    assert( vec_ );
  }

  //! iterator pointing to begin of array 
  DofIteratorType begin() {
    return DofIteratorType(*this, 0);
  }
  
  //! const iterator pointing to begin of array 
  ConstDofIteratorType begin() const {    
    return ConstDofIteratorType(*this, 0);
  }
  
  //! iterator pointing to end of array 
  DofIteratorType end() {
    return DofIteratorType(*this, size_);
  }
  
  //! const iterator pointing to end of array 
  ConstDofIteratorType end() const {    
    return ConstDofIteratorType(*this, size_);
  }

  //! return number of enties of array 
  int size () const { return size_; }  

  //! return reference to entry i
  T& operator [] ( int i )       
  { 
    assert( ((i<0) || (i>=size ()) ? (std::cout << std::endl << i << " i|size " << size() << std::endl, 0) : 1));
    return vec_[i]; 
  }
  
  //! return reference to const entry i
  const T& operator [] ( int i ) const 
  { 
    assert( ((i<0) || (i>=size()) ? (std::cout << std::endl << i << " i|size " << size() << std::endl, 0) : 1));
    return vec_[i]; 
  }

  //! assign arrays 
  ThisType& operator= (const ThisType & org)
  {
    assert(org.size_ >= size() );
    std::memcpy(vec_, org.vec_, size_ * sizeof(T));
    return *this;
  }
 
  //! operator +=  
  ThisType& operator += (const ThisType & org)
  {
    assert(org.size_ >= size() );
    const int s = size();
    const T * ov = org.vec_;
    for(int i=0; i<s; ++i) vec_[i] += ov[i];
    return *this;
  }
 
  //! operator -=  
  ThisType& operator -= (const ThisType& org)
  {
    assert(org.size() >= size() );
    const int s = size();
    const T * ov = org.vec_;
    for(int i=0; i<s; ++i) vec_[i] -= ov[i];
    return *this;
  }
 
  //! axpy operation  
  void axpy (const ThisType& org, const T scalar)
  {
    const int s = size();
    const T * ov = org.vec_;
    for(int i=0; i<s; ++i) vec_[i] += scalar*ov[i];
  }
 
  //! set all entries to zero 
  void clear () 
  {
    const int s = size();
    for(int i=0; i<s; ++i) vec_[i] = 0;
  }
 
  //! move memory from old to new destination 
  void memmove(const int length, const int oldStartIdx, const int newStartIdx) 
  {
    void * dest = ((void *) (&vec_[newStartIdx]));
    const void * src = ((const void *) (&vec_[oldStartIdx]));
    std::memmove(dest, src, length * sizeof(T));
  }
 
  //! operator = assign all entrys with value t 
  ThisType& operator= (const T t)
  {
    const int s = size();
    for(int i=0; i<s; ++i) vec_[i] = t;
    return *this;
  }

  //! Comparison operator
  //! The comparison operator checks for object identity, i.e. if this and
  //! other are the same objects in memory rather than containing the same data
  bool operator==(const ThisType& other) const 
  {
    return vec_ == other.vec_;
  }

  //! return leak pointer for usage in BLAS routines 
  T* leakPointer() { return vec_; }
  //! return leak pointer for usage in BLAS routines 
  const T* leakPointer() const { return vec_; }

  //! read and write xdr 
  bool processXdr(XDR *xdrs)
  {
    if(xdrs != 0)
    {
      int len = size_;
      xdr_int( xdrs, &len );
      assert(size_ <= len);

      xdr_vector(xdrs,(char *) vec_,size_, sizeof(T) ,(xdrproc_t)xdr_double);
      return true;
    }
    else
      return false;
  } 

  //! print array 
  void print(std::ostream& s) const 
  {
    s << "Print StaticArray(addr = "<< this << ") (size = " << size_ << ")\n";
    for(int i=0; i<size(); ++i)
    {
      s << vec_[i] << "\n";
    }
  }
};

// specialisations of axpy 
template <>
inline void StaticArray<double>::axpy(const ThisType& org, const double scalar)
{
#if HAVE_BLAS
  DuneCBlas :: daxpy( size() , scalar, org.vec_, 1 , vec_, 1);
#else 
  const int s = size();
  const double * ov = org.vec_;
  for(int i=0; i<s; ++i) vec_[i] += scalar * ov[i];
#endif
}
 
// specialisations of clear 
template <>
inline void StaticArray<int>::clear()
{
  std::memset(vec_, 0 , size() * sizeof(int));
}
template <>
inline void StaticArray<double>::clear()
{
  std::memset(vec_, 0 , size() * sizeof(double));
}
 
//! specialisation for int 
template <>
inline bool StaticArray<int>::processXdr(XDR *xdrs)
{
  typedef int T;
  if(xdrs != 0)
  {
    int len = size_;
    xdr_int( xdrs, &len );
    assert(size_ <= len);
    xdr_vector(xdrs,(char *) vec_,size_, sizeof(T) ,(xdrproc_t)xdr_int);
    return true;
  }
  else
    return false;
}

//! specialisation for double 
template <>
inline bool StaticArray<double>::processXdr(XDR *xdrs)
{
  typedef double T;
  
  if(xdrs != 0)
  {
    int len = size_;
    xdr_int( xdrs, &len );
    assert( (size_ > len) ? (std::cout << size_ << " s|l " << len << "\n" ,0 ): 1);

    xdr_vector(xdrs,(char *) vec_,size_, sizeof(T) ,(xdrproc_t)xdr_double);
    return true;
  }
  else
    return false;
}

/*! 
 MutableArray is the array that a discrete functions sees. If a discrete
 function is created, then it is signed in by the function space and the
 return value is a MemObject. This MemObject contains a MutableArrayMemory
 which is then as reference given to the MutableArray of the DiscreteFunction. 
 The MutableArray is only a wrapper class for MutableArrayMemory where we dont know
 the type of the dofs only the size of one dof. 
 Therefore we have this wrapper class for cast to the right type.
*/
template <class T, class AllocatorType>
class MutableArray : public StaticArray<T>
{
private:
  typedef MutableArray<T, AllocatorType> ThisType;
  typedef StaticArray<T> BaseType;

  // actual capacity of array
  int memSize_;
 
  // make new memory memFactor larger 
  double memoryFactor_;

  MutableArray(const MutableArray&);
public:
  //! create array of length size
  //! if size is <= 0 then vec of with memory 1 and size 0 is created (parallel runs)
  MutableArray(const int size = 0) 
    : BaseType((size < 0) ? 0 : size, 
               AllocatorType :: malloc ((size<=0) ? 1 : size))
    , memSize_((size<=0) ? 1 : size) 
    , memoryFactor_(1.0)
  {
  }
  
  //! set memory factor
  void setMemoryFactor(const double memFactor)
  {
    memoryFactor_ = memFactor;
  }

  //! Destructor 
  ~MutableArray() 
  {
    if( this->vec_ ) AllocatorType :: free ( this->vec_ );
  }

  //! return number of total enties of array 
  int capacity () const { return memSize_; }  
 
  //! resize vector with new size nsize
  //! if nsize is smaller then actual memSize, size is just set to new value
  void resize ( int nsize )
  {
    assert(nsize >= 0);

    // just set size and do not change memory in this case 
    if( (nsize <= memSize_) && (nsize> (memSize_/2)) ) 
    {
      this->size_ = nsize;
      return ;
    }

    // nsize is the minimum needed size of the vector 
    // we double this size to reserve some memory and minimize
    // reallocations 
    assert( memoryFactor_ >= 1.0 );
    const double overEstimate = memoryFactor_ * nsize;
    const int nMemSize = (int) overEstimate;
    assert( nMemSize >= nsize );

    // reallocate memory 
    this->vec_ = AllocatorType :: realloc (this->vec_,memSize_,nMemSize);

    this->size_ = nsize;
    memSize_ = nMemSize;
  }

  //! return size of vector in bytes 
  int usedMemorySize() const 
  {
    return memSize_ * sizeof(T);
  } 
};

} // end namespace Dune 
#endif
