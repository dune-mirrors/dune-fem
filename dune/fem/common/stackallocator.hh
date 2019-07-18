#ifndef DUNE_FEM_COMMON_STACKALLOCATOR_HH
#define DUNE_FEM_COMMON_STACKALLOCATOR_HH

#include <cassert>
#include <new>
#include <stack>
#include <utility>

namespace Dune
{

  namespace Fem
  {

    // UninitializedObjectStack
    // ------------------------

    struct UninitializedObjectStack
      : public std::stack< void * >
    {
      explicit UninitializedObjectStack ( std::size_t objectSize )
        : objectSize_( objectSize )
      {}

      UninitializedObjectStack ( const UninitializedObjectStack &other )
        : std::stack< void * >(),
          objectSize_( other.objectSize_ )
      {}

      UninitializedObjectStack &operator= ( const UninitializedObjectStack &other )
      {
        if( objectSize_ != other.objectSize_ )
          clear();
        objectSize_ = other.objectSize_;
        return *this;
      }

      ~UninitializedObjectStack () { clear(); }

      void clear ()
      {
        for( ; !empty(); pop() )
          ::operator delete( top() );
      }

      std::size_t objectSize () const { return objectSize_; }

      void resize ( std::size_t newSize ) { clear(); objectSize_ = newSize; }

    private:
      std::size_t objectSize_;
    };



    // StackAllocator
    // --------------

    template< class T, class S = UninitializedObjectStack * >
    struct StackAllocator
    {
      typedef T value_type;

      typedef T *pointer;
      typedef const T*const_pointer;

      typedef T &reference;
      typedef const T &const_reference;

      typedef std::size_t size_type;
      typedef std::ptrdiff_t difference_type;

      template< class U >
      struct rebind { typedef StackAllocator< U, S > other; };

      typedef UninitializedObjectStack Stack;
      typedef S StackPtr;

      explicit StackAllocator ( StackPtr stack ) : stack_( stack ) {}

      template< class U >
      StackAllocator ( const StackAllocator< U, S > &other ) : stack_( other.stack_ ) {}

      template< class U >
      StackAllocator ( StackAllocator< U, S > &&other ) : stack_( std::move( other.stack_ ) ) {}

      StackAllocator ( const StackAllocator &other ) : stack_( other.stack_ ) {}
      StackAllocator ( StackAllocator && other ) : stack_( other.stack_ ) {}

      pointer address ( reference x ) const { return &x; }
      const_pointer address ( const_reference x ) const { return &x; }

      pointer allocate ( size_type n, typename rebind< void >::other::const_pointer hint = nullptr )
      {
        assert( n <= max_size() );
        if( !stack().empty() )
        {
          pointer p = (pointer) stack().top();
          stack().pop();
          return p;
        }
        else
          return (pointer) ::operator new( stack().objectSize() );
      }

      void deallocate ( pointer p, size_type n )
      {
        assert( n <= max_size() );
        stack().push( p );
      }

      size_type max_size () const { return stack().objectSize() / sizeof( T ); }

      template< class... Args >
      void construct ( pointer p, Args &&... args )
      {
        assert( p );
        new( p ) T( std::forward< Args >( args )... );
      }

      void destroy ( pointer p ) { p->~T(); }

    private:
      template< class, class >
      friend struct StackAllocator;

      const Stack &stack () const { return *stack_; }
      Stack &stack () { return *stack_; }

      StackPtr stack_;
    };


    template<class S>
    struct StackAllocator<void, S>
    {
      typedef void value_type;

      typedef void *pointer;
      typedef const void*const_pointer;

      typedef std::size_t size_type;
      typedef std::ptrdiff_t difference_type;

      template< class U >
      struct rebind { typedef StackAllocator< U, S > other; };

      typedef UninitializedObjectStack Stack;
    };

  } // namespace

} //namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_STACKALLOCATOR_HH
