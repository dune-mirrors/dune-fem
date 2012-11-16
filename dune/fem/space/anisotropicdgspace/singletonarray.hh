#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SINGLETONARRAY_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SINGLETONARRAY_HH

// C++ includes
#include <cassert>
#include <cstddef>
#include <utility>

// dune-common includes
#include <dune/common/array.hh>
#include <dune/common/nullptr.hh>

// dune-fem includes
#include <dune/fem/misc/threadmanager.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me.
*/


namespace Dune
{

  namespace Fem
  {

    // SingletonArrayDefaultFactory
    // ----------------------------

    template< class Key, class Object >
    struct SingletonArrayDefaultFactory;

    template< class Object >
    struct SingletonArrayDefaultFactory< std::size_t, Object >
    {
      static Object *createObject ( const std::size_t &i ) { return new Object( i ); }

      static void deleteObject ( Object *object ) { delete object; }

      static std::size_t index ( const std::size_t &i ) { return i; }
    };

    template< class Key, class Object >
    struct SingletonArrayDefaultFactory 
    {
      static Object *createObject ( const Key &key ) { return new Object( key ); }

      static void deleteObject( Object *object ) { delete object; }

      static std::size_t index ( const Key &key ) { return std::size_t( key ); }
    };



    // SingletonArray
    // --------------

    template< class T, int N, class Key = std::size_t, 
              class Factory = SingletonArrayDefaultFactory< std::size_t, T > >
    class SingletonArray
    {
      typedef SingletonArray< T, N, Factory > ThisType;

      typedef Factory FactoryType;

    public:
      typedef std::size_t size_type;
      typedef T value_type;
      typedef Key key_type;

    private:
      SingletonArray ()
      {
        for( size_type i = 0; i < N; ++i )
          entries_[ i ] = std::pair< value_type *, size_type>( nullptr, 0 );
      }

      SingletonArray ( const ThisType & );

    public:
      static SingletonArray &instance ()
      {
        static ThisType instance;
        return instance;
      }

      inline static value_type &get ( const key_type &key )
      {
        assert( Fem::ThreadManager::singleThreadMode() );

        // get entry
        const size_type index = FactoryType::index( key );
        std::pair< value_type *, size_type > &entry = instance().entries_[ index ];

        // increment reference counter, if entry already exists
        if( entry.first )
        {
          assert( entry.second > 0 );
          ++entry.second;
          return *(entry.first);
        }
        
        // otherwise, create value
        value_type *value = FactoryType::createObject( index );
        assert( value );

        // store value
        entry = std::pair< value_type *, size_type >( value, 1 );
        return *(entry.first);
      }

      inline static void remove ( const key_type &key )
      {
        assert( Fem::ThreadManager::singleThreadMode() );
        
        // get entry
        const size_type index = FactoryType::index( key );
        std::pair< value_type *, size_type > &entry = instance().entries_[ index ];

        // check, if entry exists
        if( !entry.first )
          return;

        // get reference counter
        size_type &refCount = entry.second;
        assert( refCount > 0 );
      
        // decrement reference counter
        --refCount;
        if( refCount == 0 )
        {
          FactoryType::deleteObject( entry.first );
          entry = std::pair< value_type *, size_type>( nullptr, 0 );
        }
      }

      inline static size_type size () { return N; };

      inline static size_type max_size () { return N; };

      inline static bool empty ()
      {
        for( size_type i = 0; i < N; ++i )
        {
          std::pair< value_type *, size_type > &entry = instance().entries_[ i ];
          if( entry.first )
          {
            assert( entry.second > 0 );
            return false;
          }
        }
        return true;
      }

    private:
      array< std::pair< value_type *, size_type >, N > entries_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SINGLETONARRAY_HH
