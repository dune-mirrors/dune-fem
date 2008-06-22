#ifndef DUNE_ARRAYS_HH
#define DUNE_ARRAYS_HH

//- System includes 
#include <cassert>
#include <cstring>
#include <cstdlib>

#include <string>

//- Dune includes 
#include <dune/common/genericiterator.hh>
#include <dune/common/interfaces.hh>
#include <dune/common/exceptions.hh>

#if HAVE_BLAS 
// include BLAS for daxpy operation 
#include <dune/fem/solver/oemsolver/cblas.h>
#endif

// include xdr wrapper 
#include <dune/fem/io/file/xdrio.hh>

namespace Dune {

// forward declarations 
template <class T>
class DefaultDofAllocator;

template <class T, class AllocatorType = DefaultDofAllocator<T> >
class MutableArray;

template<class ArrayType>
struct SpecialArrayFeatures;


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
  typedef T FieldType;
  //! definition conforming to STL  
  typedef T value_type;
  
  //! definition conforming to ISTL  
  typedef T block_type;
  
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
 
  //! operator *= multiplies array with a scalar  
  ThisType& operator *= (const T scalar)
  {
    const int s = size();
    for(int i=0; i<s; ++i) vec_[i] *= scalar;
    return *this;
  }
  
  //! operator /= divides array with a scalar  
  ThisType& operator /= (const T scalar)
  {
    const T scalar_1 = (((T) 1)/scalar); 
    const int s = size();
    for(int i=0; i<s; ++i) vec_[i] *= scalar_1;
    return *this;
  }
  
  //! operator = assign all entrys to given scalar value  
  ThisType& operator= (const T scalar)
  {
    const int s = size();
    for(int i=0; i<s; ++i) vec_[i] = scalar;
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
  bool processXdr(XDRStream& xdr)
  {
    int len = size_;
    xdr.inout( len );

    // when read check size 
    if( size_ != len )
    {
      DUNE_THROW(InvalidStateException,"StaticArray::processXdr: internal size " << size_ << " and size to read " << len << " not equal!");
    }
    return processXdrVector(xdr);
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
  
protected:  
  //! read and write vector as bytes using xdr method xdr_bytes  
  bool processXdrVector(XDRStream& xdr)
  {
    // write/read all entries 
    int ret = 1;
    for(int i=0; i<size_; ++i)
    {
      ret |= xdr.inout( vec_[i] );
    }
    return (ret == 1) ? true : false;
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
  const double* ov = org.vec_;
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
  //! create array of length 0 
  MutableArray() 
    : BaseType(0,0)
    , memSize_(0) 
    , memoryFactor_(1.0)
  {
  }
  
  //! create array of length size
  MutableArray(const int size) 
    : BaseType((size < 0) ? 0 : size, 
               // only alloc memory if size > 0
               (size <= 0) ? 0 : AllocatorType :: malloc (size))
    , memSize_((size<=0) ? 0 : size) 
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
    freeMemory();
  }

  //! return number of total enties of array 
  int capacity () const { return memSize_; }  
 
  //! resize vector with new size nsize
  //! if nsize is smaller then actual memSize, size is just set to new value
  void resize ( int nsize )
  {
    // just set size if nsize is smaller than memSize but larger the
    // half of memSize 
    if( (nsize <= memSize_) && (nsize > (memSize_/2))) 
    {
      this->size_ = nsize;
      return ;
    }

    // if nsize == 0 freeMemory
    if( nsize <= 0 ) 
    {
      freeMemory();
      return ;
    }

    // reserve memory + overestimate 
    reserve( nsize );
    // set new size 
    this->size_ = nsize;
  }

  //! reserve vector size with new mSize 
  //! if mSize is smaller then actual memSize, 
  //! then nothing is done 
  void reserve ( int mSize )
  {
    // check whether we already have the mem size 
    // and if just do nothing 
    if( mSize <= memSize_ ) 
    {
      return ;
    }

    assert( memoryFactor_ >= 1.0 );
    const double overEstimate = memoryFactor_ * mSize;
    const int nMemSize = (int) overEstimate;
    assert( nMemSize >= mSize );

    if( !this->vec_ )
    {
      // allocate new memory 
      this->vec_ = AllocatorType :: malloc(nMemSize);
    }
    else 
    {
      assert( nMemSize > 0 );
      // nsize is the minimum needed size of the vector 
      // we double this size to reserve some memory and minimize
      // reallocations 
      assert( this->vec_ );

      // reallocate memory 
      this->vec_ = AllocatorType :: realloc (this->vec_,memSize_,nMemSize);
    }
    // set new mem size 
    memSize_ = nMemSize;
  }

  //! return size of vector in bytes 
  int usedMemorySize() const 
  {
    return memSize_ * sizeof(T) + sizeof(ThisType);
  } 
  
  //! read and write xdr, during read resize is done 
  //! if sizes do not match  
  bool processXdr(XDRStream& xdr)
  {
    int len = this->size();
    xdr.inout( len );

    // if actual size is smaller then resize vector  
    if( len > this->size() ) resize ( len );

    // write array 
    return this->processXdrVector(xdr);
  }

private: 
  // free memory and reset sizes 
  void freeMemory() 
  {
    if( this->vec_ ) 
    {
      AllocatorType :: free ( this->vec_ );
      this->vec_ = 0;
    }
    this->size_ = 0;
    memSize_ = 0;
  }
};

/** \brief Specialization of SpecialArrayFeatures for MutableArray */
template<class ValueType>
struct SpecialArrayFeatures<MutableArray<ValueType> >
{
  typedef MutableArray<ValueType> ArrayType;
  static size_t used(const ArrayType & array)  
  {
    return array.usedMemorySize();
  }
  static void setMemoryFactor(ArrayType & array, const double memFactor) 
  {
    array.setMemoryFactor(memFactor);
  }

  static void memMoveBackward(ArrayType& array, const int length,
            const int oldStartIdx, const int newStartIdx)
  {
    //array.memmove(length,oldStartIdx,newStartIdx);
    // get new end of block which is offSet + (length of block - 1) 
    int newIdx = newStartIdx + length - 1; 
    // copy all entries backwards 
    for(int oldIdx = oldStartIdx+length-1; oldIdx >= oldStartIdx; --oldIdx, --newIdx )
    {
      // move value to new location 
      array[newIdx] = array[oldIdx];
#ifndef NDEBUG
      // for debugging purpose 
      array[oldIdx ] = 0.0;
#endif
    }
  }
  static void memMoveForward(ArrayType& array, const int length,
            const int oldStartIdx, const int newStartIdx)
  {
    //array.memmove(length,oldStartIdx,newStartIdx);
    const int upperBound = oldStartIdx + length;
    // get new off set that should be smaller then old one
    int newIdx = newStartIdx;
    for(int oldIdx = oldStartIdx; oldIdx<upperBound; ++oldIdx, ++newIdx )
    {
      // copy to new location 
      array[newIdx] = array[oldIdx];
#ifndef NDEBUG 
      // for debugging issues only 
      array[oldIdx] = 0.0;
#endif
    }
  }
};

} // end namespace Dune 
#endif
