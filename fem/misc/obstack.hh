#ifndef DUNE_FEM_OBJECTSTACK_HH
#define DUNE_FEM_OBJECTSTACK_HH

#include <obstack.h>
#include <iostream>

#include <dune/fem/misc/arrayallocator.hh>


namespace Dune
{

  #define obstack_chunk_alloc &std :: malloc
  #define obstack_chunk_free &std :: free
  #define obstack_alloc_failed_handler = &ObStack :: allocFailed;

  class ObStack
  {
  private:
    typedef ObStack ThisType;

  protected:
    struct obstack info_;

  public:
    inline ObStack ()
    {
      obstack_init( &info_ );
      // align on 8 byte boundaries
      obstack_alignment_mask( &info_ ) = 7;
    }

    inline ObStack ( unsigned int chunkSize )
    {
      obstack_init( &info_ );
      obstack_chunk_size( &info_ ) = chunkSize;
      obstack_alignment_mask( &info_ ) = 7;
    }

    inline ~ObStack ()
    {
      // free the entire object stack
      obstack_free( &info_, 0 );
    }
    
    inline void *allocate ( unsigned int size )
    {
      return obstack_alloc( &info_, size );
    }

    inline void free ( void *ptr )
    {
      if( ptr != 0 )
        obstack_free( &info_, ptr );
    }

  private:
    static inline void allocFailed ()
    {
      DUNE_THROW( OutOfMemoryError, "ObStack went out of memory." );
    }
  };

  #undef obstack_chunk_alloc
  #undef obstack_chunk_free



  template< class ElementImp >
  class ObStackArrayAllocator;



  template< class ElementImp >
  struct ObStackArrayAllocatorTraits
  {
    typedef ElementImp ElementType;

    typedef ElementImp *ElementPtrType;

    typedef ObStackArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class ElementImp >
  class ObStackArrayAllocator
  : public ArrayAllocatorDefault< ObStackArrayAllocatorTraits< ElementImp > >
  {
  public:
    typedef ElementImp ElementType;

    typedef ObStackArrayAllocatorTraits< ElementType > TraitsType;

  private:
    typedef ObStackArrayAllocator< ElementType > ThisType;
    typedef ArrayAllocatorDefault< TraitsType > BaseType;

  public:
    typedef typename TraitsType :: ElementPtrType ElementPtrType;

  public:
    inline void allocate ( unsigned int size,
                           ElementPtrType &array ) const
    {
      array = (ElementPtrType)obStack().allocate( size * sizeof( ElementType ) );
      assert( array != 0 );
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      obStack().free( array );
      array = 0;
    }

    inline void reallocate ( unsigned int oldSize,
                             unsigned int newSize,
                             ElementPtrType &array ) const
    {
      DUNE_THROW( NotImplemented,
                  "Reallocating an object on an object stack is not possible." );
    }

  private:
    static inline ObStack& obStack ()
    {
      static ObStack stack( 2 << 18 );
      return stack;
    }
  };

}

#endif
