#ifndef DUNE_FEM_SINGLETON_HH
#define DUNE_FEM_SINGLETON_HH

//- System includes
#include <cassert>
#include <algorithm>
#include <memory>
#include <typeindex>
#include <unordered_map>
#include <vector>

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
      protected:
        // item to be stored in storage list
        struct Item {  virtual ~Item() {} };

        typedef std::shared_ptr< Item > WeakPointerType;
        typedef std::unique_ptr< Item > PointerType;
        typedef std::type_index         KeyType;

        typedef std::pair< std::unordered_map< KeyType, std::shared_ptr< Item > >, std::vector<PointerType> >  StorageType;

        // delete for singletons deleting each singleton object
        // in reverse order of creation
        struct SingletonDeleter
        {
          void operator()(StorageType* storage) const
          {
            // delete singletons in reverse order
            std::for_each(storage->second.rbegin(), storage->second.rend(),
                          [](PointerType& item) { item.reset(); });

            storage->second.clear();
            storage->first.clear();
          }
        };

      public:
        typedef std::unique_ptr<StorageType, SingletonDeleter> StoragePointer;

      private:
        static StoragePointer storage_;

      protected:
        static StorageType& getStorage()
        {
          if(! storage_ )
          {
            storage_.reset( new StorageType() );
          }

          return *storage_;
        }
      };
    } // end namespace detail

    /** \brief return singleton instance of given Object type.
     */
    template< class Object >
    class Singleton : public detail::SingletonStorage
    {
      typedef detail::SingletonStorage BaseType;
      typedef typename BaseType::StorageType      StorageType;
      typedef typename BaseType::Item             Item;
      typedef typename BaseType::PointerType      PointerType;
      typedef typename BaseType::WeakPointerType  WeakPointerType;

      using BaseType::getStorage;

      // item to be created containing the correct object
      struct ItemWrapper : public Item
      {
        Object obj_;
        template <class... Args>
        ItemWrapper(Args &&... args) : obj_(std::forward< Args >( args )...)
        {}
      };
      typedef ItemWrapper ItemWrapperType;

      // null delete for shared_ptr hash map
      struct NullDeleter
      {
        void operator()(Item *p) const {}
      };

    public:
      /** \brief return singleton instance of given Object type.
       *  \param args Possible constructor arguments for object when created for first time
       */
      template <class... Args>
      static Object& instance(Args &&... args)
      {
        // capture reference as static variable to avoid map search later on
        static Object& inst = getObject(std::forward< Args >( args )...);
        return inst;
      }

    protected:
      // placing variables as static inside functions only works with gcc
      static const bool placeStaticVariableInline = false ;

      /** \brief return singleton instance of given Object type.
       */
      template <class... Args>
      static Object& getObject(Args &&... args)
      {
        // this way of creating static variables only works with gcc, not with clang
        if constexpr ( placeStaticVariableInline )
        {
          static Object obj( std::forward< Args >( args )...);
          return obj;
        }
        else
        {
          // get storage reference (see base class)
          StorageType& storage = getStorage();

          // get pointer of singleton objects belonging to hash id
          auto& ptr = storage.first[ std::type_index(typeid(Object)) ];
          // if pointer has not been set, create object and set pointer
          if( ! ptr )
          {
            assert( Fem::ThreadManager::singleThreadMode() );

            // create object in vector for later correct deletion order
            storage.second.emplace_back( PointerType(new ItemWrapperType(std::forward< Args >( args )...) ) );

            // create pointer to object in hash map for later use
            ptr = WeakPointerType( storage.second.back().operator->(), NullDeleter() );
          }

          // return object reference
          assert( dynamic_cast< ItemWrapperType* > (ptr.operator->()) );
          return static_cast< ItemWrapperType& > (*ptr).obj_;
        }
      }

    };

  } // namespace Fem

} // namespace Dune

#endif //  #ifndef DUNE_FEM_SINGLETONLIST_HH
