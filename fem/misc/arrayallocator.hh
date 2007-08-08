#ifndef DUNE_FEM_ARRAYALLOCATOR_HH
#define DUNE_FEM_ARRAYALLOCATOR_HH

#include <cstdlib>

namespace Dune
{
 
  template< class TraitsImp >
  class ArrayAllocatorInterface
  {
  public:
    typedef TraitsImp TraitsType;

  private:
    typedef ArrayAllocatorInterface< TraitsType > ThisType;

  public:
    typedef ThisType ArrayAllocatorInterfaceType;
    
    typedef typename TraitsType :: ArrayAllocatorType ArrayAllocatorType;
    
    typedef typename TraitsType :: ElementType ElementType;

    typedef typename TraitsType :: ElementPtrType ElementPtrType;

  public:
    inline ArrayAllocatorInterface ()
    {
      typedef CompileTimeChecker< Conversion< ArrayAllocatorType, ThisType > :: exists >
        __Array_Allocator_Implementation_Must_Be_Derived_From_Interface__;
    }
    
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

  protected:
    inline const ArrayAllocatorType &asImp () const
    {
      return static_cast< const ArrayAllocatorType& >( *this );
    }

    inline ArrayAllocatorType &asImp ()
    {
      return static_cast< ArrayAllocatorType& >( *this );
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

}

#endif
