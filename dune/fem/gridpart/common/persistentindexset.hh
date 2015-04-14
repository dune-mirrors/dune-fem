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



    // PersistentIndexSet
    // ------------------

    /** \brief please doc me
     *
     *  \tparam  Grid  grid type
     *  \tparam  Implementation  index set type
     */
    template< class Traits, template< class > class Base >
    class PersistentIndexSet
      : public Base< Traits >,
        public PersistentIndexSetInterface
    {
      typedef Base< Traits > BaseType;

    protected:
      using BaseType::impl;

      /** \brief grid type */
      typedef typename Traits::GridType GridType;
      /** \brief dof manager type */
      typedef DofManager< GridType > DofManagerType;

      explicit PersistentIndexSet ( const GridType &grid )
        : grid_( grid ),
          dofManager_( DofManagerType::instance( grid ) ),
          counter_( 0 )
      {
        dofManager_.addIndexSet( impl() );
      }

    public:
      using BaseType::read;
      using BaseType::write;

      ~PersistentIndexSet ()
      {
        dofManager_.removeIndexSet( impl() );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::backup */
      void backup () const
      {
        if( needsBackupRestore() )
          write( PersistenceManager::backupStream() );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::backup */
      void restore ()
      {
        if( needsBackupRestore() )
          read( PersistenceManager::restoreStream() );
      }

      /** \copydoc Dune::Fem::PersistentIndexSetInterface::addBackupRestore */
      void addBackupRestore () override final { ++counter_; }

      /** \copydoc Dune::Fem::PersistentIndexSetInterface::removeBackupRestore */
      void removeBackupRestore () override final { --counter_; }

    private:
      bool needsBackupRestore () const { return counter_ > 0; }

    protected:
      const GridType &grid_;
      DofManagerType &dofManager_;

    private:
      int counter_ ;
    };



    // PersistentConsecutiveIndexSet
    // -----------------------------

    template< class Traits >
    class PersistentConsecutiveIndexSet
      : public PersistentIndexSet< Traits, ConsecutiveIndexSet >
    {
      typedef PersistentIndexSet< Traits, ConsecutiveIndexSet > BaseType;

    protected:
      explicit PersistentConsecutiveIndexSet ( const typename BaseType::GridType &grid )
        : BaseType( grid )
      {}
    };



    // PersistentAdaptiveIndexSet
    // --------------------------

    template< class Traits >
    class PersistentAdaptiveIndexSet
      : public PersistentIndexSet< Traits, AdaptiveIndexSet >
    {
      typedef PersistentIndexSet< Traits, AdaptiveIndexSet > BaseType;

    protected:
      explicit PersistentAdaptiveIndexSet ( const typename BaseType::GridType &grid )
        : BaseType( grid )
      {}
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_TEST_PERSISTENTINDEXSET_HH
