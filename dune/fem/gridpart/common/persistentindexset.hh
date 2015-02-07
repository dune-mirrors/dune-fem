#ifndef DUNE_FEM_GRIDPART_TEST_PERSISTENTINDEXSET_HH
#define DUNE_FEM_GRIDPART_TEST_PERSISTENTINDEXSET_HH

#include <type_traits>
#include <utility>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/io/file/persistencemanager.hh>
#include <dune/fem/io/streams/streams.hh>
#include <dune/fem/space/common/dofmanager.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class Traits >
    class PersistentIndexSet;
    template< class Traits >
    class PersistentConsecutiveIndexSet;
    template< class Traits >
    class PersistentAdaptiveIndexSet;



    // PersistentIndexSetInterface
    // ---------------------------

    /** \brief virtual base class for persistent index sets
     *
     */
    struct PersistentIndexSetInterface
    {
      virtual ~PersistentIndexSetInterface () {}

      /** \brief please doc me */
      virtual void addBackupRestore () = 0;

      /** \brief please doc me */
      virtual void removeBackupRestore () = 0;
    };



    namespace Capabilities
    {

      // isPersistentIndexSet
      // --------------------

#ifndef DOXYGEN

      template< class IndexSet, bool value = std::is_base_of< PersistentIndexSetInterface, IndexSet >::type::value >
      struct __isPersistentIndexSet;

      template< class IndexSet >
      struct __isPersistentIndexSet< IndexSet, true >
      {
        static const bool v = true;

        static constexpr PersistentIndexSetInterface* map ( IndexSet &indexSet )
        {
          return static_cast< PersistentIndexSetInterface * >( &indexSet );
        }
      };

      template< class IndexSet >
      struct __isPersistentIndexSet< IndexSet, false >
      {
        static const bool v = false;

        static constexpr PersistentIndexSetInterface* map ( IndexSet & ) noexcept
        {
          return nullptr;
        }
      };

#endif // #ifndef DOXYGEN

      /** \brief capability for persistent index sets
       *
       *  \tparam  IndexSet  an index set type
       *
       *  \note default value is \b true if IndexSet is derived from
       *        PersistentIndexSetInterface
       */
      template< class IndexSet >
      struct isPersistentIndexSet
        : public __isPersistentIndexSet< IndexSet >
      {
      private:
        typedef __isPersistentIndexSet< IndexSet > BaseType;

      public:
        /** \brief please doc me */
        static const bool v = BaseType::v;

        /** \brief please doc me */
        static constexpr PersistentIndexSetInterface* map ( IndexSet &indexSet ) noexcept
        {
          return BaseType::map( indexSet );
        }
      };

#ifndef DOXYGEN

      template< class IndexSet >
      struct isPersistentIndexSet< const IndexSet >
        : public isPersistentIndexSet< IndexSet >
      {};

#endif // #ifndef DOXYGEN

    } // namespace Capabilities



    // PersistentIndexSetBase
    // ----------------------

    /** \brief please doc me
     *
     *  \tparam  Grid  grid type
     *  \tparam  Implementation  index set type
     */
    template< class Grid, class Implementation >
    class PersistentIndexSetBase
      : public PersistentObject
    {
    protected:
      typedef Grid GridType;
      typedef DofManager< GridType > DofManagerType;

      explicit PersistentIndexSetBase ( const GridType &grid )
        : grid_( grid ),
          dofManager_( DofManagerType::instance( grid ) ),
          counter_( 0 )
      {
        dofManager_.addIndexSet( impl() );
      }

    public:
      /** \brief please doc me */
      void backup () const override final
      {
        if( needsBackupRestore() )
          impl().write( PersistenceManager::backupStream() );
      }

      /** \brief please doc me */
      void restore () override final
      {
        if( needsBackupRestore() )
          impl().read( PersistenceManager::backupStream() );
      }

      /** \brief please doc me */
      void addBackupRestore () override final { ++counter_; }

      /** \brief please doc me */
      void removeBackupRestore () override final { --counter_; }

      /** \brief please doc me */
      template< class T >
      void write ( OutStreamInterface< T > &stream )
      {
        impl().write( stream );
      }

      /** \brief please doc me */
      template< class T >
      void read ( InStreamInterface< T > &stream )
      {
        impl().read( stream );
      }

    private:
      bool needsBackupRestore () const { return counter_ > 0; }

      Implementation &impl ()
      {
        return static_cast< Implementation & >( *this );
      }

      const Implementation &impl () const
      {
        return static_cast< const Implementation & >( *this );
      }

    protected:
      const GridType &grid_;
      DofManagerType &dofManager_;

    public:
      int counter_ ;
    };



    // PersistentConsecutiveIndexSet
    // -----------------------------

    template< class Traits >
    class PersistentConsecutiveIndexSet
      : public ConsecutiveIndexSet< Traits >,
        public PersistentIndexSetBase< typename Traits::GridType, typename Traits::IndexSetType >
    {
      typedef PersistentIndexSetBase< typename Traits::GridType, typename Traits::IndexSetType > BaseType;

    protected:
      using BaseType::BaseType;
    };



    // PersistentAdaptiveIndexSet
    // --------------------------

    template< class Traits >
    class PersistentAdaptiveIndexSet
      : public AdaptiveIndexSet< Traits >,
        public PersistentIndexSetBase< typename Traits::GridType, typename Traits::IndexSetType >
    {
      typedef PersistentIndexSetBase< typename Traits::GridType, typename Traits::IndexSetType > BaseType;

    protected:
      using BaseType::BaseType;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_PERSISTENTINDEXSET_HH
