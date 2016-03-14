#ifndef DUNE_FEM_ARRAYALLOCATOR_HH
#define DUNE_FEM_ARRAYALLOCATOR_HH

#include <cstdlib>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  namespace Fem
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

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ARRAYALLOCATOR_HH
