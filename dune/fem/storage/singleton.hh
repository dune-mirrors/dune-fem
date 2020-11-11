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

      public:
        template <class Object>
        struct ItemWrapper : public Item
        {
          Object obj_;
          template <class... Args>
          ItemWrapper(Args &&... args) : obj_(std::forward< Args >( args )...)
          {}
        };

        struct NullDeleter
        {
          void operator()(Item *p) {}
        };

        typedef std::unique_ptr< Item > PointerType;
        typedef std::type_index         KeyType;

        typedef std::unordered_map< KeyType, PointerType > StorageType;

      private:
        static StorageType* storage_;

        // placing variables as static inside functions only works with gcc
        static const bool placeStaticVariableInline = false ;

      public:
        /** \brief return singleton instance of given Object type.
         */
        template <class Object, class... Args>
        static Object& instance(Args &&... args)
        {
          // this way of creating static variables only works with gcc, not with clang
          if constexpr ( placeStaticVariableInline )
          {
            static Object obj( std::forward< Args >( args )...);
            return obj;
          }
          else
          {
            if(! storage_ )
              storage_ = new StorageType();

            StorageType& storage = *storage_;
            //std::cout << "Accessing Object " << typeid(Object).name() << std::endl;
            //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;

            typedef ItemWrapper< Object > ItemWrapperType;
            PointerType& ptr = storage[ std::type_index(typeid(Object)) ];
            if( ! ptr )
            {
              assert( Fem::ThreadManager::singleThreadMode() );
              ptr.reset( new ItemWrapperType(std::forward< Args >( args )...) );
              //std::cout << "Create Object " << typeid(Object).name() << std::endl;
              //std::cout << "typeindex = " << std::type_index(typeid(Object)).hash_code() << std::endl;
            }
            assert( dynamic_cast< ItemWrapperType* > (ptr.operator->()) );
            return static_cast< ItemWrapperType& > (*ptr).obj_;
          }
        }
      };
    }

    /** \brief return singleton instance of given Object type.
     */
    template< class Object >
    struct DUNE_EXPORT Singleton
    {
      /** \brief return singleton instance of given Object type.
       *  \param args Possible constructor arguments for object when created for first time
       */
      template <class... Args>
      static Object& instance(Args &&... args)
      {
        // catch reference as static variable to avoid map search later on
        static Object& inst = detail::SingletonStorage::template instance< Object, Args... > (std::forward< Args >( args )...);
        return inst;
      }
    };

  } // namespace Fem

} // namespace Dune

#endif //  #ifndef DUNE_FEM_SINGLETONLIST_HH
