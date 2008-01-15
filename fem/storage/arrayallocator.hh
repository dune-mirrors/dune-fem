#ifndef DUNE_FEM_ARRAYALLOCATOR_HH
#define DUNE_FEM_ARRAYALLOCATOR_HH

#include <cstdlib>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{
 
  template< class TraitsImp >
  class ArrayAllocatorInterface
  : public BartonNackmanInterface< ArrayAllocatorInterface< TraitsImp >,
                                   typename TraitsImp :: ArrayAllocatorType >
  {
  public:
    typedef TraitsImp TraitsType;

    typedef typename TraitsType :: ArrayAllocatorType ArrayAllocatorType;

  private:
    typedef ArrayAllocatorInterface< TraitsType > ThisType;
    typedef BartonNackmanInterface< ThisType, ArrayAllocatorType > BaseType;

  public:
    typedef ThisType ArrayAllocatorInterfaceType;
    
    typedef typename TraitsType :: ElementType ElementType;

    typedef typename TraitsType :: ElementPtrType ElementPtrType;

  protected:
    using BaseType :: asImp;
    
  public:
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().allocate( size, array ) );
    }

    inline void free ( ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().free( array ) );
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().reallocate( oldSize, newSize, array ) );
    }

    inline void reserve ( unsigned int newSize,
                          ElementPtrType &array ) const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION
        ( asImp().reserve( newSize, array ) );
    }
  };



  template< class TraitsImp >
  class ArrayAllocatorDefault
  : public ArrayAllocatorInterface< TraitsImp >
  {
  public:
    typedef TraitsImp TraitsType;
    
  private:
    typedef ArrayAllocatorDefault< TraitsType > ThisType;
    typedef ArrayAllocatorInterface< TraitsType > BaseType;

  public:
    using BaseType :: allocate;
    using BaseType :: free;

  public:
    typedef typename TraitsType :: ElementType ElementType;

    typedef typename TraitsType :: ElementPtrType ElementPtrType;

  public:
    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      ElementPtrType newArray;
      allocate( newSize, newArray );
      
      const unsigned int copySize = std :: min( oldSize, newSize );
      for( unsigned int i = 0; i < copySize; ++i )
        newArray[ i ] = array[ i ];

      free( array );
      array = newArray;
    }

    // by default, do nothing
    inline void reserve ( unsigned int newSize,
                          ElementPtrType &array ) const
    {
    }
  };



  template< class ElementImp >
  class StandardArrayAllocator;

  template< class ElementImp >
  class CArrayAllocator;



  // Choose a default array allocator
  #ifndef USE_CARRAYALLOCATOR
    template< class ElementType >
    class DefaultArrayAllocator
    : public StandardArrayAllocator< ElementType >
    {
    };
  #else
    template< class ElementType >
    class DefaultArrayAllocator
    : public CArrayAllocator< ElementType >
    {
    };
  #endif


  
  template< class ElementImp >
  struct StandardArrayAllocatorTraits
  {
    typedef ElementImp ElementType;

    typedef ElementImp *ElementPtrType;

    typedef StandardArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class ElementImp >
  class StandardArrayAllocator
  : public ArrayAllocatorDefault< StandardArrayAllocatorTraits< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef StandardArrayAllocatorTraits< ElementType > TraitsType;

  private:
    typedef StandardArrayAllocator< ElementType > ThisType;
    typedef ArrayAllocatorDefault< TraitsType > BaseType;

  public:
    typedef typename TraitsType :: ElementPtrType ElementPtrType;
    
  public:
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      if( size > 0 )
      {
        array = new ElementType[ size ];
        assert( array != 0 );
      }
      else
        array = 0;
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      if( array != 0 )
      {
        delete[]( array );
        array = 0;
      }
    }
  };



  template< class ElementImp >
  struct CArrayAllocatorTraits
  {
    typedef ElementImp ElementType;

    typedef ElementImp *ElementPtrType;

    typedef CArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class ElementImp >
  class CArrayAllocator
  : public ArrayAllocatorDefault< CArrayAllocatorTraits< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef CArrayAllocatorTraits< ElementType > TraitsType;

  private:
    typedef CArrayAllocator< ElementType > ThisType;
    typedef ArrayAllocatorDefault< TraitsType > BaseType;

  public:
    typedef typename TraitsType :: ElementPtrType ElementPtrType;
    
  public:
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      if( size > 0 )
      {
        array = (ElementPtrType)malloc( size * sizeof( ElementType ) );
        assert( array != 0 );
      } else
        array = 0;
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      if( array != 0 )
      {
        std :: free( array );
        array = 0;
      }
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      if( newSize == 0 )
      {
        std :: free( array );
        array = 0;
      }
      else
        array = (ElementPtrType)realloc( array, newSize * sizeof( ElementType ) );
    }
  };



  template< class ElementImp, template< class > class WrappedArrayAllocatorImp >
  class ArrayOverAllocator;



  template< class ElementImp, template< class > class WrappedArrayAllocatorImp >
  class ArrayOverAllocatorElementPointer
  {
  public:
    typedef ElementImp ElementType;

    typedef WrappedArrayAllocatorImp< ElementType > WrappedArrayAllocatorType;

    typedef typename WrappedArrayAllocatorType :: ElementPtrType ElementPtrType;

  private:
    typedef ArrayOverAllocatorElementPointer
      < ElementType, WrappedArrayAllocatorImp >
      ThisType;

    friend class ArrayOverAllocator< ElementType, WrappedArrayAllocatorImp >;

  protected:
    ElementPtrType ptr_;
    unsigned int size_;

  public:
    inline ArrayOverAllocatorElementPointer ()
    : ptr_( 0 ),
      size_( 0 )
    {
    }

    inline ArrayOverAllocatorElementPointer ( const ElementPtrType ptr,
                                              const unsigned int size )
    : ptr_( ptr ),
      size_( size )
    {
    }

    inline ArrayOverAllocatorElementPointer ( const ThisType &other )
    : ptr_( other.ptr_ ),
      size_( other.size_ )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      ptr_ = other.ptr_;
      size_ = other.size_;
    }

#if 0
    inline operator const ElementPtrType () const
    {
      return ptr_;
    }
#endif

    inline ElementType &operator* () const
    {
      return *ptr_;
    }

    inline ElementType &operator[] ( const unsigned int index ) const
    {
      assert( index < size_ );
      return ptr_[ index ];
    }
  };



  template< class ElementImp, template< class > class WrappedArrayAllocatorImp >
  struct ArrayOverAllocatorTraits
  {
    typedef ElementImp ElementType;

    typedef WrappedArrayAllocatorImp< ElementType > WrappedArrayAllocatorType;

    typedef ArrayOverAllocatorElementPointer
      < ElementType, WrappedArrayAllocatorImp >
      ElementPtrType;

    typedef ArrayOverAllocator< ElementType, WrappedArrayAllocatorImp >
      ArrayAllocatorType;
  };



  template< class ElementImp, template< class > class WrappedArrayAllocatorImp >
  class ArrayOverAllocator
  : public ArrayAllocatorDefault
    < ArrayOverAllocatorTraits< ElementImp, WrappedArrayAllocatorImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef ArrayOverAllocatorTraits< ElementType, WrappedArrayAllocatorImp >
      TraitsType;

  private:
    typedef ArrayOverAllocator< ElementType, WrappedArrayAllocatorImp > ThisType;
    typedef ArrayAllocatorDefault< TraitsType > BaseType;

  public:
    typedef typename TraitsType :: ElementPtrType ElementPtrType;

    typedef typename TraitsType :: WrappedArrayAllocatorType
      WrappedArrayAllocatorType;

  protected:
    WrappedArrayAllocatorType allocator_;
    unsigned int memFactor_;

  public:
    inline explicit ArrayOverAllocator ( const unsigned int memFactor = 1152 )
    : allocator_(),
      memFactor_( memFactor )
    {
    }

    inline explicit ArrayOverAllocator ( const double memFactor )
    : allocator_(),
      memFactor_( (int)(memFactor * 1024) )
    {
    }


    inline explicit
    ArrayOverAllocator ( const WrappedArrayAllocatorType &allocator,
                         const unsigned int memFactor = 1152 )
    : allocator_( allocator ),
      memFactor_( memFactor )
    {
    }

    inline ArrayOverAllocator ( const WrappedArrayAllocatorType &allocator,
                                const double memFactor )
    : allocator_( allocator ),
      memFactor_( (int)(memFactor * 1024) )
    {
    }

    inline ArrayOverAllocator ( const ThisType &other )
    : allocator_( other.allocator_ ),
      memFactor_( other.memFactor_ )
    {
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      allocator_ = other.allocator_;
      memFactor_ = other.memFactor_;
    }

    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      array.size_ = (size * memFactor_) / 1024;
      allocator_.allocate( array.size_, array.ptr_ );
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      allocator_.free( array.ptr_ );
      array.size_ = 0;
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      const unsigned int newAllocSize = (newSize * memFactor_) / 1024;
      if( array.size_ < newSize )
        reserve( newAllocSize );
    }

    inline void reserve ( unsigned int newSize,
                          ElementPtrType &array ) const
    {
      if( newSize > array.size_ )
      {
        allocator_.reallocate( array.size_, newSize, array.ptr_ );
        array.size_ = newSize;
      }
    }
  };



  template< class ElementType >
  class DefaultArrayOverAllocator
  : public ArrayOverAllocator< ElementType, DefaultArrayAllocator >
  {
  };

}

#endif
