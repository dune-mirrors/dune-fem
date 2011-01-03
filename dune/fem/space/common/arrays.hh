#ifndef DUNE_ARRAYS_HH
#define DUNE_ARRAYS_HH

//- System includes 
#include <cassert>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

#include <string>

//- Dune includes 
#include <dune/common/genericiterator.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>

#if HAVE_BLAS 
// include BLAS for daxpy operation 
#include <dune/fem/solver/oemsolver/cblas.h>
#endif
#if DUNE_VERSION_NEWER( DUNE_GRID, 2, 1, 0 )
#include <dune/grid/alugrid/common/interfaces.hh>
#else
#include <dune/grid/alugrid/interfaces.hh>
#endif

// include xdr wrapper 
#include <dune/fem/io/streams/streams.hh>

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
    T* p = new T [ nmemb ] ;
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
    assert(oldMem);
    assert(nmemb > 0);
    T* p = malloc(nmemb);
    const size_t copySize = std::min( oldSize, nmemb );
    std::copy( oldMem, oldMem+copySize, p );
    free (oldMem);
    return p;
  }
};

//! allocator for simple structures like int, double and float
//! using the C malloc,free, and realloc 
template <typename T> 
struct SimpleDofAllocator 
{
  //! allocate array of nmemb objects of type T
  static T* malloc (size_t nmemb)
  {
    assert(nmemb > 0);
    T* p = (T *) std::malloc(nmemb * sizeof(T));
    assert(p);
    return p;
  }

  //! release memory previously allocated with malloc member
  static void free (T* p)
  {
    assert(p);
    std::free(p);
  }
  
  //! allocate array of nmemb objects of type T
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
struct DefaultDofAllocator<double> : public SimpleDofAllocator< double > 
{
};

template <>
struct DefaultDofAllocator< float > : public SimpleDofAllocator< float > 
{
};

template <>
struct DefaultDofAllocator< int > : public SimpleDofAllocator< int > 
{
};

template <>
struct DefaultDofAllocator< size_t > : public SimpleDofAllocator< size_t > 
{
};

template <>
struct DefaultDofAllocator< char > : public SimpleDofAllocator< char > 
{
};

template <>
struct DefaultDofAllocator< bool > : public SimpleDofAllocator< bool > 
{
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

  // pointer to mem
  T * vec_;

  // size of array 
  size_t size_;

  StaticArray(const StaticArray&);
public:
  typedef T FieldType;
  //! definition conforming to STL  
  typedef T value_type;
  
  //! definition conforming to ISTL  
  typedef T block_type;
  
  //! DofIterator
  typedef GenericIterator<ThisType, T> DofIteratorType;
  
  //! make compatible with std::vector 
  typedef DofIteratorType iterator ;

  //! Const DofIterator
  typedef GenericIterator<const ThisType, const T> ConstDofIteratorType;

  //! make compatible with std::vector 
  typedef ConstDofIteratorType const_iterator ;

  //! create array of length size and store vec as pointer to memory 
  explicit StaticArray(const size_t size, T* vec) 
    : vec_(vec) 
    , size_(size)
  {
    assert( size_ >= 0 );
  }

  //! create array of length size and store vec as pointer to memory 
  explicit StaticArray(const size_t size, const T* vec) 
    : vec_( const_cast< T * > (vec) ) 
    , size_(size)
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
  size_t size () const { return size_; }  

private:
  void assertIndex ( const size_t i ) const
  {
#ifndef NDEBUG
    if( i >= size() )
    {
      std::cerr << std::endl;
      std::cerr << "Error in StaticArray: Index out of Range: " << i << std::endl;
      std::cerr << "                      Size of array: " << size() << std::endl;
      abort();
    }
#endif
  }

public:
  //! return reference to entry i
  T& operator [] ( const size_t i )       
  { 
    assertIndex( i );
    return vec_[i]; 
  }
  
  //! return reference to const entry i
  const T& operator [] ( const size_t i ) const 
  {
    assertIndex( i );
    return vec_[i]; 
  }

  //! assign arrays 
  ThisType& operator= (const ThisType & org)
  {
    assert(org.size_ >= size() );
    assert( ( size_ > 0 ) ? vec_ != 0 : true );
    // copy the entries 
    std::copy(org.vec_, org.vec_ + size_, vec_ );
    return *this;
  }
 
  //! operator +=  
  ThisType& operator += (const ThisType & org)
  {
    assert(org.size_ >= size() );
    const size_t s = size();
    const T * ov = org.vec_;
    for(size_t i=0; i<s; ++i) vec_[i] += ov[i];
    return *this;
  }
 
  //! operator -=  
  ThisType& operator -= (const ThisType& org)
  {
    assert(org.size() >= size() );
    const size_t s = size();
    const T * ov = org.vec_;
    for(size_t i=0; i<s; ++i) vec_[i] -= ov[i];
    return *this;
  }
 
  //! operator *= multiplies array with a scalar  
  ThisType& operator *= (const T scalar)
  {
    const size_t s = size();
    for(size_t i=0; i<s; ++i) vec_[i] *= scalar;
    return *this;
  }
  
  //! operator /= divides array with a scalar  
  ThisType& operator /= (const T scalar)
  {
    const T scalar_1 = (((T) 1)/scalar); 
    const size_t s = size();
    for(size_t i=0; i<s; ++i) vec_[i] *= scalar_1;
    return *this;
  }
  
  //! operator = assign all entrys to given scalar value  
  ThisType& operator= (const T scalar)
  {
    const size_t s = size();
    for(size_t i=0; i<s; ++i) vec_[i] = scalar;
    return *this;
  }

  //! axpy operation  
  void axpy (const ThisType& org, const T scalar)
  {
    const size_t s = size();
    const T * ov = org.vec_;
    for(size_t i=0; i<s; ++i) vec_[i] += scalar*ov[i];
  }
 
  //! set all entries to zero 
  void clear () 
  {
    const size_t s = size();
    for(size_t i=0; i<s; ++i) vec_[i] = 0;
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

  //! write to  stream 
  template <class StreamTraits> 
  bool write(OutStreamInterface< StreamTraits >& out) const 
  {
    out << size_;
    for(size_t i=0; i<size_; ++i)
    {
      out << vec_[i];
    }
    return true;
  }

  //! write to  stream 
  template <class StreamTraits> 
  bool read(InStreamInterface< StreamTraits >& in) 
  {
    size_t len; 
    in >> len;
    // when read check size 
    if( size_ != len )
    {
      DUNE_THROW(InvalidStateException,"StaticArray::read: internal size " << size_ << " and size to read " << len << " not equal!");
    }

    for(size_t i=0; i<size_; ++i)
    {
      in >> vec_[i];
    }
    return true;
  }

  //! print array 
  void print(std::ostream& s) const 
  {
    s << "Print StaticArray(addr = "<< this << ") (size = " << size_ << ")\n";
    for(size_t i=0; i<size(); ++i)
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
  const size_t s = size();
  const double* ov = org.vec_;
  for(size_t i=0; i<s; ++i) vec_[i] += scalar * ov[i];
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
protected:
  typedef MutableArray<T, AllocatorType> ThisType;
  typedef StaticArray<T> BaseType;

  using BaseType :: size_ ;
  using BaseType :: vec_ ;

  // make new memory memFactor larger 
  double memoryFactor_;

  // actual capacity of array
  size_t memSize_;
 
public:
  using BaseType :: size ;

  //! create array of length 0 
  MutableArray() 
    : BaseType(0, (T *) 0)
    , memoryFactor_(1.0)
    , memSize_(0) 
  {
  }

  //! copy constructor
  MutableArray(const MutableArray& other) 
    : BaseType(0, (T *) 0),
      memoryFactor_(1.0),
      memSize_(0)
  {
    // assign vector 
    *this = other;
  }

  //! create array of length size
  MutableArray(const size_t size) 
    : BaseType(size, 
               // only alloc memory if size > 0
               ((T *) (size == 0) ? 0 : AllocatorType :: malloc (size)))
    , memoryFactor_(1.0)
    , memSize_(size) 
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
  size_t capacity () const { return memSize_; }  
 
  //! assign arrays 
  ThisType& operator= (const ThisType & org)
  {
    resize( org.size_ );
    memoryFactor_ = org.memoryFactor_;
    assert( ( size_ > 0 ) ? vec_ != 0 : true );
    std::copy(org.vec_, org.vec_ + size_, vec_ );
    return *this;
  }
 
  //! resize vector with new size nsize
  //! if nsize is smaller then actual memSize, size is just set to new value
  void resize ( size_t nsize )
  {
    // just set size if nsize is smaller than memSize but larger the
    // half of memSize 
    if( (nsize <= memSize_) && (nsize >= (memSize_/2)) ) 
    {
      size_ = nsize;
      return ;
    }

    // if nsize == 0 freeMemory
    if( nsize == 0 ) 
    {
      freeMemory();
      return ;
    }

    // reserve or shrink to memory + overestimate 
    adjustMemory( nsize );
    // set new size 
    size_ = nsize;
  }

  //! reserve vector size with new mSize 
  //! if mSize is smaller then actual memSize, 
  //! then nothing is done 
  void reserve ( size_t mSize )
  {
    // check whether we already have the mem size 
    // and if just do nothing 
    if( mSize <= memSize_ ) 
    {
      return ;
    }

    // adjust memory accordingly 
    adjustMemory( mSize );
  }

  //! return size of vector in bytes 
  size_t usedMemorySize() const 
  {
    return memSize_ * sizeof(T) + sizeof(ThisType);
  } 

protected: 
  //! adjust the memory 
  void adjustMemory( size_t mSize )
  {
    assert( memoryFactor_ >= 1.0 );
    const double overEstimate = memoryFactor_ * mSize;
    const size_t nMemSize = (size_t) std::ceil( overEstimate );
    assert( nMemSize >= mSize );

    if( !vec_ )
    {
      // allocate new memory 
      vec_ = AllocatorType :: malloc(nMemSize);
    }
    else 
    {
      assert( nMemSize > 0 );
      // nsize is the minimum needed size of the vector 
      // we double this size to reserve some memory and minimize
      // reallocations 
      assert( vec_ );

      // reallocate memory 
      vec_ = AllocatorType :: realloc (vec_,memSize_,nMemSize);
    }

    // set new mem size 
    memSize_ = nMemSize;
  }

  // free memory and reset sizes 
  void freeMemory() 
  {
    if( vec_ ) 
    {
      AllocatorType :: free ( vec_ );
      vec_ = 0;
    }
    size_ = 0;
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

  static void memMoveBackward(ArrayType& array, 
                              const size_t length,
                              const size_t oldStartIdx, 
                              const size_t newStartIdx)
  {
    assert( newStartIdx >= oldStartIdx );
    //array.memmove(length,oldStartIdx,newStartIdx);
    // get new end of block which is offSet + (length of block - 1) 
    size_t newIdx = newStartIdx + length - 1; 
    assert( newIdx < array.size() );
    // copy all entries backwards 
    for(size_t oldIdx = oldStartIdx + length-1; oldIdx >= oldStartIdx; --oldIdx, --newIdx )
    {
      assert( oldIdx < array.size() );
      // move value to new location 
      array[newIdx] = array[oldIdx];
#ifndef NDEBUG
      // for debugging purpose 
      array[oldIdx ] = 0.0;
#endif
    }
  }
  static void memMoveForward(ArrayType& array, 
                             const size_t length,
                             const size_t oldStartIdx, 
                             const size_t newStartIdx)
  {
    assert( newStartIdx <= oldStartIdx );
    //array.memmove(length,oldStartIdx,newStartIdx);
    const size_t upperBound = oldStartIdx + length;
    // get new off set that should be smaller then old one
    size_t newIdx = newStartIdx;
    for(size_t oldIdx = oldStartIdx; oldIdx<upperBound; ++oldIdx, ++newIdx )
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
