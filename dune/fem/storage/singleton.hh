#ifndef DUNE_FEM_SINGLETON_HH
#define DUNE_FEM_SINGLETON_HH

//- System includes
#include <cassert>
#include <iostream>
#include <memory>
#include <typeindex>
#include <utility>
#include <unordered_map>

#include <dune/common/visibility.hh>

//- dune-fem includes
#include <dune/fem/misc/threads/threadmanager.hh>

namespace Dune
{
  namespace Fem
  {
    namespace detail
    {
      class SingletonStorage
      {
      public:
        typedef std::shared_ptr< void > PointerType;
        typedef std::type_index          KeyType;

        typedef std::unordered_map< KeyType, PointerType > StorageType;

      private:
        static StorageType storage_;

        template <class Object>
        struct Deleter
        {
          void operator()(void* p) const
          {
            delete (Object *) p;
          }
        };

      public:
        /** \brief return singleton instance of given Object type.
         */
        template <class Object>
        DUNE_EXPORT
        static Object& instance()
        {
          PointerType& ptr = storage_[ std::type_index(typeid(Object)) ];
          if( ! ptr )
          {
            assert( Fem::ThreadManager::singleThreadMode() );
            //std::cout << "Create Object " << typeid(Object).name() << std::endl;
            //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;
            Object* obj = new Object();
            ptr = PointerType( (void *) obj, Deleter< Object > () );
          }
          return *((Object *) ptr.operator->());
        }
      };
    }

    template< class Object >
    struct Singleton
    {
      /** \brief return singleton instance of given Object type.
       */
      DUNE_EXPORT
      static Object& instance()
      {
        // TODO: introduce cmake variable
#if 1
        // def DUNE_FEM_USE_HASED_BASED_SINGLETON
        return detail::SingletonStorage::template instance< Object > ();
#else
        // traditional way of creating static variable
        // that does not work with older clang compilers
        static Object obj;
        return obj;
#endif
      }
    };

  } // namespace Fem

} // namespace Dune

#endif //  #ifndef DUNE_FEM_SINGLETONLIST_HH
