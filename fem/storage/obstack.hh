#ifndef DUNE_FEM_OBJECTSTACK_HH
#define DUNE_FEM_OBJECTSTACK_HH

#include <obstack.h>
#include <iostream>

#include <dune/fem/storage/arrayallocator.hh>


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
      // align on 16 byte boundaries
      obstack_alignment_mask( &info_ ) = 15;
    }

    inline ObStack ( unsigned int chunkSize )
    {
      obstack_init( &info_ );
      obstack_chunk_size( &info_ ) = chunkSize;
      // align on 16 byte boundaries
      obstack_alignment_mask( &info_ ) = 15;
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



  template< class Element >
  class ObStackArrayAllocator;



  template< class Element >
  struct ObStackArrayAllocatorTraits
  {
    typedef Element ElementType;

    typedef Element *ElementPtrType;

    typedef ObStackArrayAllocator< ElementType > ArrayAllocatorType;
  };



  template< class Element >
  class ObStackArrayAllocator
  : public ArrayAllocatorDefault< ObStackArrayAllocatorTraits< Element > >
  {
  public:
    typedef Element ElementType;

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
      if( size > 0 )
      {
        array = (ElementPtrType)obStack().allocate( size * sizeof( ElementType ) );
        assert( array != 0 );
      }
      else
        array = 0;
    }
  
    inline void free ( ElementPtrType &array ) const
    {
      if( array != 0 )
      {
        obStack().free( array );
        array = 0;
      }
    }

  private:
    static inline ObStack& obStack ()
    {
      static ObStack stack( 2 << 20 );
      return stack;
    }
  };

}

#endif
