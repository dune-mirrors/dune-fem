#ifndef DUNE_FEM_SINGLETON_HH
#define DUNE_FEM_SINGLETON_HH

//- System includes
#include <cassert>
#include <iostream>
#include <memory>
#include <typeindex>
#include <utility>
#include <map>
#include <unordered_map>
#include <vector>

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
        struct Item
        {
          virtual ~Item() {}
        };

        template <class Object>
        struct ItemWrapper : public Item
        {
          Object obj_;
        };

      public:
        struct NullDeleter
        {
          void operator()(Item *p) {}
        };

        //typedef std::unique_ptr< Item, NullDeleter > PointerType;
        typedef std::unique_ptr< Item > PointerType;
        typedef std::type_index         KeyType;

        //typedef std::map< KeyType, PointerType > StorageType;
        //typedef std::vector< PointerType > StorageType;
        typedef std::unordered_map< KeyType, PointerType > StorageType;

      private:
        static StorageType* storage_;

      public:
        /** \brief return singleton instance of given Object type.
         */
        template <class Object>
        DUNE_EXPORT
        static Object& instance()
        {
          if(! storage_ )
            storage_ = new StorageType();

          StorageType& storage = *storage_;
#if 0
          static Object obj;
          return obj;
#else
          //std::cout << "Accessing Object " << typeid(Object).name() << std::endl;
          //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;

          typedef ItemWrapper< Object > ItemWrapperType;
          PointerType& ptr = storage[ std::type_index(typeid(Object)) ];
          if( ! ptr )
          {
            assert( Fem::ThreadManager::singleThreadMode() );
            ptr.reset( new ItemWrapperType() );
            //std::cout << "Create Object " << typeid(Object).name() << std::endl;
            //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;
          }
          assert( dynamic_cast< ItemWrapperType* > (ptr.operator->()) );
          return static_cast< ItemWrapperType& > (*ptr).obj_;
#endif
        }
      };
    }

    /** \brief return singleton instance of given Object type.
     */
    template< class Object >
    struct DUNE_EXPORT Singleton
    {
      //static std::unique_ptr< Object > instance_;

      /** \brief return singleton instance of given Object type.
       */
      DUNE_EXPORT
      static Object& instance()
      {
        static Object& inst = detail::SingletonStorage::template instance< Object > ();
        return inst;
        //if( ! instance_ )
        //  instance_.reset( new Object() );
        //return *instance_;
        // forward to non-templated class
        //return detail::SingletonStorage::template instance< Object > ();
      }
    };

    //template <class Object> std::unique_ptr< Object > DUNE_EXPORT Singleton< Object >::instance_;

  } // namespace Fem

} // namespace Dune

#endif //  #ifndef DUNE_FEM_SINGLETONLIST_HH
