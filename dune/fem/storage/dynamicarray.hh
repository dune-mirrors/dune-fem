#ifndef DUNE_FEM_DYNAMICARRAY_HH
#define DUNE_FEM_DYNAMICARRAY_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <type_traits>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>
#include <dune/fem/common/utility.hh>

namespace Dune
{

  namespace Fem
  {
  // forward declaration
  template< class K > class StaticArray;

  //! oriented to the STL Allocator funtionality
  template <typename T>
  class StandardArrayAllocator
    : public std::allocator< T >
  {
    typedef std::allocator< T > BaseType;
  public:
#if __cplusplus <= 201703L
    typedef typename BaseType :: pointer    pointer ;
#else
    typedef T* pointer;
#endif
    typedef typename BaseType :: size_type  size_type;

    pointer allocate( size_type n )
    {
      return new T[ n ];
    }

    void deallocate( pointer p, size_type n )
    {
      delete [] p;
    }

    pointer reallocate ( pointer oldMem, size_type oldSize, size_type n )
    {
      assert(oldMem);
      pointer p = allocate( n );
      const size_type copySize = std::min( oldSize, n );
      std::copy( oldMem, oldMem+copySize, p );
      deallocate( oldMem, oldSize );
      return p;
    }
  };

  //! allocator for simple structures like int, double and float
  //! using the C malloc, free and realloc
  template <typename T>
  class PODArrayAllocator : public std::allocator< T >
  {
    static_assert( Std::is_pod< T > :: value, "T is not POD" );
    typedef std::allocator< T > BaseType;
  public:
    PODArrayAllocator() = default;

#if __cplusplus <= 201703L
    typedef typename BaseType :: pointer    pointer ;
#else
    typedef T* pointer;
#endif
    typedef typename BaseType :: size_type  size_type;
    typedef typename BaseType :: value_type value_type;

    //! allocate array of nmemb objects of type T
    pointer allocate( size_type n )
    {
      pointer p = static_cast< pointer > (std::malloc(n * sizeof(value_type)));
      assert(p);
      return p;
    }

    //! release memory previously allocated with malloc member
    void deallocate( pointer p, size_type n )
    {
      assert(p);
      std::free(p);
    }

    //! allocate array of nmemb objects of type T
    pointer reallocate (pointer oldMem, size_type oldSize , size_type n)
    {
      assert(oldMem);
      pointer p = static_cast< pointer > (std::realloc(oldMem , n*sizeof(value_type)));
      assert(p);
      return p;
    }
  };

  template <class T, class AllocatorType = typename std::conditional< Std::is_pod< T > :: value,
                                             PODArrayAllocator< T >,
                                             StandardArrayAllocator< T > > :: type >
  class DynamicArray;

  } // end namespace Fem

  // specialization of DenseMatVecTraits for StaticArray
  template< class K >
  struct DenseMatVecTraits< Fem::StaticArray< K > >
  {
    typedef Fem::StaticArray< K > derived_type;
    typedef K* container_type;
    typedef K value_type;
    typedef std::size_t size_type;
  };

  template< class K >
  struct FieldTraits< Fem::StaticArray< K > >
  {
    typedef typename FieldTraits< K >::field_type field_type;
    typedef typename FieldTraits< K >::real_type real_type;
  };

  // specialization of DenseMatVecTraits for DynamicArray
  template< class K >
  struct DenseMatVecTraits< Fem::DynamicArray< K > > : public DenseMatVecTraits< Fem::StaticArray< K > >
  {
  };

  template< class K >
  struct FieldTraits< Fem::DynamicArray< K > > : public FieldTraits< Fem::StaticArray< K > >
  {
  };

  namespace Fem
  {


  /** \brief An implementation of DenseVector which uses a C-array of fixed size as storage
    *
    * \tparam T is the field type (use float, double, complex, etc)
    */
  template <class T>
  class StaticArray : public DenseVector< StaticArray< T > >
  {
    typedef StaticArray< T > ThisType;
    typedef DenseVector< ThisType > BaseType;

  public:
    typedef typename BaseType::size_type size_type;

    typedef typename BaseType::value_type value_type;
    typedef value_type FieldType;

    typedef typename DenseMatVecTraits< ThisType >::container_type DofStorageType;

    StaticArray(const ThisType&) = delete;

    //! create array of length size and store vec as pointer to memory
    explicit StaticArray(size_type size, const value_type* vec)
      : vec_( const_cast< DofStorageType > (vec) )
      , size_(size)
    {
    }

    //! return size of array
    size_type size () const
    {
      return size_;
    }

    //! random access operator
    value_type& operator [] ( size_type i )
    {
      assert( i < size_ );
      return vec_[i];
    }

    //! random access operator
    const value_type& operator [] ( size_type i ) const
    {
      assert( i < size_ );
      return vec_[i];
    }

    //! copy assignament
    ThisType& operator= (const ThisType & org)
    {
      assert(org.size_ >= size() );
      assert( ( size_ > 0 ) ? vec_ != nullptr : true );
      std::copy(org.vec_, org.vec_ + size_, vec_ );
      return *this;
    }

    //! set all entries to 0
    void clear ()
    {
      std::fill( vec_, vec_+size_, value_type(0) );
    }

    //! move memory from old to new destination
    void memmove(size_type length, size_type oldStartIdx, size_type newStartIdx)
    {
      void * dest = ((void *) (&vec_[newStartIdx]));
      const void * src = ((const void *) (&vec_[oldStartIdx]));
      std::memmove(dest, src, length * sizeof(value_type));
    }

    //! comparison operator: checks for object identity, i.e. if this and
    //! other are the same objects in memory rather than containing the same data
    bool operator==(const ThisType& other) const
    {
      return vec_ == other.vec_;
    }

    //! return pointer to data
    value_type* data()
    {
      return vec_;
    }

    //! return pointer to data
    const value_type* data() const
    {
      return vec_;
    }

  protected:
    DofStorageType vec_;
    size_type size_;
  };


  /** \brief An implementation of DenseVector which uses a C-array of dynamic size as storage
    *
    * \tparam T is the field type (use float, double, complex, etc)
    * \tparam Allocator is the allocator type
    */
  template <class T, class Allocator>
  class DynamicArray : public StaticArray<T>
  {
  public:
    typedef Allocator  AllocatorType;
  protected:
    typedef DynamicArray<T, AllocatorType> ThisType;
    typedef StaticArray<T> BaseType;

    using BaseType :: size_ ;
    using BaseType :: vec_ ;

  public:
    using BaseType :: size ;

    typedef typename BaseType::size_type size_type;
    typedef typename BaseType::value_type value_type;

    //! copy constructor
    DynamicArray(const ThisType& other)
      : BaseType(0, nullptr)
      , memoryFactor_(1.0)
      , memSize_(0)
      , allocator_( other.allocator_ )
    {
      assign( other );
    }

    //! create array of length size with initialized values
    explicit DynamicArray(size_type size,
                          const value_type& value,
                          AllocatorType allocator = AllocatorType() )
      : BaseType(size, (size == 0) ? nullptr : allocator.allocate(size) )
      , memoryFactor_(1.0)
      , memSize_(size)
      , allocator_( allocator )
    {
      if( size_ > 0 )
      {
        std::fill( vec_, vec_+size_, value );
      }
    }

    //! create array of length size without initializing the values
    explicit DynamicArray(size_type size = 0,
                          AllocatorType allocator = AllocatorType() )
      : BaseType(size, (size == 0) ? nullptr : allocator.allocate(size) )
      , memoryFactor_(1.0)
      , memSize_(size)
      , allocator_( allocator )
    {
    }

    //! set memory factor
    void setMemoryFactor(double memFactor)
    {
      memoryFactor_ = memFactor;
      assert( memoryFactor_ >= 1.0 );
    }

    //! destructor
    ~DynamicArray()
    {
      freeMemory();
    }

    //! return number of total enties of array
    size_type capacity () const
    {
      return memSize_;
    }

    //! assign arrays
    void assign (const ThisType & org)
    {
      memoryFactor_ = org.memoryFactor_;
      assert( memoryFactor_ >= 1.0 );

      resize( org.size_ );
      assert( ( size_ > 0 ) ? vec_ != nullptr : true );
      std::copy(org.vec_, org.vec_ + size_, vec_ );
    }

    //! assign arrays
    ThisType& operator= (const ThisType & org)
    {
      assign( org );
      return *this;
    }

    //! resize vector with new size nsize
    //! if nsize is smaller then actual memSize, size is just set to new value
    void resize ( size_type nsize )
    {
      // only initialize value if we are not using a POD type
      doResize( nsize, ! Std::is_pod< value_type >::value );
    }

    //! resize vector with new size nsize
    //! if nsize is smaller then actual memSize, size is just set to new value
    void resize ( size_type nsize, const value_type& value )
    {
      doResize( nsize, true, value );
    }

    void doResize( size_type nsize, bool initializeNewValues, const value_type& value = value_type() )
    {
      // just set size if nsize is smaller than memSize but larger the half of memSize
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
      adjustMemory( nsize, initializeNewValues, value );
      // set new size
      size_ = nsize;
    }

    //! reserve vector size with new mSizeif
    //! if mSize is smaller then actual memSize, then nothing is done
    void reserve ( size_type mSize )
    {
      // check whether we already have the mem size and if just do nothing
      if( mSize <= memSize_ )
      {
        return ;
      }

      // adjust memory accordingly
      adjustMemory( mSize, false );
    }

    //! return size of vector in bytes
    size_type usedMemorySize() const
    {
      return memSize_ * sizeof(value_type) + sizeof(ThisType);
    }

  protected:
    //! adjust the memory
    void adjustMemory( size_type mSize, bool initializeNewValues, const value_type& value = value_type() )
    {
      assert( memoryFactor_ >= 1.0 );
      const double overEstimate = memoryFactor_ * mSize;
      const size_type nMemSize = (size_type) std::ceil( overEstimate );
      assert( nMemSize >= mSize );

      if( vec_ == nullptr )
      {
        // allocate new memory
        vec_ = allocator_.allocate( nMemSize );
        if( initializeNewValues )
        {
          std::fill( vec_, vec_+nMemSize, value );
        }
      }
      else
      {
        assert( nMemSize > 0 );
        assert( vec_ );

        // reallocate memory
        vec_ = allocator_.reallocate (vec_, memSize_, nMemSize);
        if( nMemSize > memSize_ && initializeNewValues )
        {
          std::fill( vec_+memSize_, vec_+nMemSize, value );
        }
      }

      // set new mem size
      memSize_ = nMemSize;
    }

    // free memory and reset sizes
    void freeMemory()
    {
      if( vec_ != nullptr )
      {
        allocator_.deallocate( vec_, memSize_ );
        vec_ = nullptr;
      }
      size_ = 0;
      memSize_ = 0;
    }

    double memoryFactor_;
    size_type memSize_;
    AllocatorType allocator_;
  };

  } // namespace Fem

} // namespace Dune
#endif // #ifndef DUNE_FEM_DYNAMICARRAY_HH
