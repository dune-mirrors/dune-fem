#ifndef DUNE_FEM_ARRAYALLOCATOR_HH
#define DUNE_FEM_ARRAYALLOCATOR_HH

#include <cstdlib>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{
 
  template< class Traits >
  class ArrayAllocatorInterface
  : public BartonNackmanInterface< ArrayAllocatorInterface< Traits >,
                                   typename Traits :: ArrayAllocatorType >
  {
    typedef ArrayAllocatorInterface< Traits > ThisType;
    typedef BartonNackmanInterface
      < ThisType, typename Traits :: ArrayAllocatorType >
      BaseType;

  public:
    typedef typename Traits :: ArrayAllocatorType ArrayAllocatorType;

    typedef ThisType ArrayAllocatorInterfaceType;
    
    typedef typename Traits :: ElementType ElementType;
    typedef typename Traits :: ElementPtrType ElementPtrType;

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



  template< class Traits >
  class ArrayAllocatorDefault
  : public ArrayAllocatorInterface< Traits >
  {
    typedef ArrayAllocatorDefault< Traits > ThisType;
    typedef ArrayAllocatorInterface< Traits > BaseType;

  public:
    using BaseType :: allocate;
    using BaseType :: free;

  public:
    typedef typename Traits :: ElementType ElementType;
    typedef typename Traits :: ElementPtrType ElementPtrType;

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
    {}
  };



  template< class Element >
  class StandardArrayAllocator;

  template< class Element >
  class CArrayAllocator;



  // Choose a default array allocator
  #ifndef USE_CARRAYALLOCATOR
    template< class Element >
    class DefaultArrayAllocator
    : public StandardArrayAllocator< Element >
    {};
  #else
    template< class Element >
    class DefaultArrayAllocator
    : public CArrayAllocator< Element >
    {};
  #endif


  
  template< class Element >
  struct StandardArrayAllocatorTraits
  {
    typedef Element ElementType;
    typedef Element *ElementPtrType;

    typedef StandardArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class Element >
  class StandardArrayAllocator
  : public ArrayAllocatorDefault< StandardArrayAllocatorTraits< Element > >
  {
    typedef StandardArrayAllocator< Element > ThisType;
    typedef StandardArrayAllocatorTraits< Element > Traits;
    typedef ArrayAllocatorDefault< Traits > BaseType;

  public:
    typedef typename Traits :: ElementType ElementType;
    typedef typename Traits :: ElementPtrType ElementPtrType;
    
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



  template< class Element >
  struct CArrayAllocatorTraits
  {
    typedef Element ElementType;

    typedef Element *ElementPtrType;

    typedef CArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class Element >
  class CArrayAllocator
  : public ArrayAllocatorDefault< CArrayAllocatorTraits< Element > >
  {
    typedef CArrayAllocator< Element > ThisType;
    typedef CArrayAllocatorTraits< Element > Traits;
    typedef ArrayAllocatorDefault< Traits > BaseType;
    
  public:
    typedef typename Traits :: Element ElementType;
    typedef typename Traits :: ElementPtrType ElementPtrType;
    
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



#if 0
  template< class Element, template< class > class WrappedArrayAllocator >
  class ArrayOverAllocator;



  template< class Element, template< class > class WrappedArrayAllocator >
  class ArrayOverAllocatorElementPointer
  {
    typedef ArrayOverAllocatorElementPointer< Element, WrappedArrayAllocator >
      ThisType;
    
    friend class ArrayOverAllocator< Element, WrappedArrayAllocator >;
    
  public:
    typedef Element ElementType;

    typedef WrappedArrayAllocator< ElementType > WrappedArrayAllocatorType;

    typedef typename WrappedArrayAllocatorType :: ElementPtrType ElementPtrType;

  protected:
    ElementPtrType ptr_;
    unsigned int size_;

  public:
    inline ArrayOverAllocatorElementPointer ()
    : ptr_( 0 ),
      size_( 0 )
    {}

    inline ArrayOverAllocatorElementPointer ( const ElementPtrType ptr,
                                              const unsigned int size )
    : ptr_( ptr ),
      size_( size )
    {}

    inline ArrayOverAllocatorElementPointer ( const ThisType &other )
    : ptr_( other.ptr_ ),
      size_( other.size_ )
    {}

    inline ThisType &operator= ( const ThisType &other )
    {
      ptr_ = other.ptr_;
      size_ = other.size_;
    }

    inline operator const ElementType * () const
    {
      return (ElementType *)ptr_;
    }

    inline operator ElementType * () const
    {
      return (ElementType *)ptr_;
    }

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



  template< class Element, template< class > class WrappedArrayAllocator >
  struct ArrayOverAllocatorTraits
  {
    typedef Element ElementType;

    typedef WrappedArrayAllocator< ElementType > WrappedArrayAllocatorType;

    typedef ArrayOverAllocatorElementPointer
      < ElementType, WrappedArrayAllocator >
      ElementPtrType;

    typedef ArrayOverAllocator< ElementType, WrappedArrayAllocator >
      ArrayAllocatorType;
  };



  template< class Element, template< class > class WrappedArrayAllocator >
  class ArrayOverAllocator
  : public ArrayAllocatorDefault
    < ArrayOverAllocatorTraits< Element, WrappedArrayAllocator > >
  {
    typedef ArrayOverAllocator< Element, WrappedArrayAllocator > ThisType;
    typedef ArrayOverAllocatorTraits< Element, WrappedArrayAllocator > Traits;
    typedef ArrayAllocatorDefault< Traits > BaseType;

  public:
    typedef typename Traits :: ElementType ElementType;

    typedef typename Traits :: ElementPtrType ElementPtrType;

    typedef typename Traits :: WrappedArrayAllocatorType
      WrappedArrayAllocatorType;

  protected:
    WrappedArrayAllocatorType allocator_;
    unsigned int memFactor_;

  public:
    inline explicit ArrayOverAllocator ( const unsigned int memFactor = 1152 )
    : allocator_(),
      memFactor_( memFactor )
    {}

    inline explicit ArrayOverAllocator ( const double memFactor )
    : allocator_(),
      memFactor_( (int)(memFactor * 1024) )
    {}


    inline explicit
    ArrayOverAllocator ( const WrappedArrayAllocatorType &allocator,
                         const unsigned int memFactor = 1152 )
    : allocator_( allocator ),
      memFactor_( memFactor )
    {}

    inline ArrayOverAllocator ( const WrappedArrayAllocatorType &allocator,
                                const double memFactor )
    : allocator_( allocator ),
      memFactor_( (int)(memFactor * 1024) )
    {}

    inline ArrayOverAllocator ( const ThisType &other )
    : allocator_( other.allocator_ ),
      memFactor_( other.memFactor_ )
    {}

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
        reserve( newAllocSize, array );
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



  template< class Element >
  class DefaultArrayOverAllocator
  : public ArrayOverAllocator< Element, DefaultArrayAllocator >
  {};
#endif

}

#endif
