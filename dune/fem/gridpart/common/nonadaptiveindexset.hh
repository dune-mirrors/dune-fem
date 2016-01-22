#ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSETWRAPPER_HH
#define DUNE_FEM_GRIDPART_COMMON_INDEXSETWRAPPER_HH

// #warning "Using Dune::Fem::NonAdaptiveIndexSet wrapper around grid part index sets"

#include <vector>

#include <dune/fem/gridpart/common/indexset.hh>

namespace Dune
{

  namespace Fem
  {

    // NonAdaptiveIndexSet
    // -------------------

    template< class IndexSet >
    class NonAdaptiveIndexSet;



    // NonAdaptiveIndexSetTraits
    // -------------------------

    template< class IndexSet >
    struct NonAdaptiveIndexSetTraits
    {
      typedef NonAdaptiveIndexSet< IndexSet > IndexSetType;

      static const int dimension = IndexSet::dimension;

      template< int codim >
      struct Codim
        : public IndexSet::template Codim< codim >
      {};

      typedef typename IndexSet::IndexType IndexType;

      typedef typename IndexSet::Types Types;
    };



    // NonAdaptiveIndexSet
    // -------------------

    /** \brief wrapper for (adaptive) index sets that disables all support for
     *         adaptivity
     */
    template< class IndexSet >
    class NonAdaptiveIndexSet
      : public AdaptiveIndexSet< NonAdaptiveIndexSetTraits< IndexSet > >
    {
      typedef AdaptiveIndexSet< NonAdaptiveIndexSetTraits< IndexSet > > BaseType;

    public:
      /** \copydoc Dune::Fem::ConsecutiveIndexSet::IndexType */
      typedef typename BaseType::IndexType IndexType;

      explicit NonAdaptiveIndexSet ( const IndexSet &indexSet )
        : indexSet_( indexSet )
      {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::types */
      typename BaseType::Types types ( int codim ) const
      {
        return indexSet().types( codim );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::contains */
      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return indexSet().contains( entity );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::size */
      IndexType size ( GeometryType type ) const
      {
        return indexSet().size( type );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::zie */
      IndexType size ( int codim ) const
      {
        return indexSet().size( codim );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::index */
      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::index */
      template< int codim >
      IndexType index ( const typename BaseType::template Codim< codim >::Entity &entity ) const
      {
        return indexSet().template index< codim >( entity );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::subIndex */
      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int cd ) const
      {
        return subIndex< Entity::codimension >( entity, i, cd );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::subIndex */
      template< int codim >
      IndexType subIndex ( const typename BaseType::template Codim< codim >::Entity &entity, int i, unsigned int cd ) const
      {
        return indexSet().template subIndex< codim >( entity, i, cd );
      }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::resize */
      static void resize () {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::compress */
      static constexpr bool compress () noexcept { return false; }

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::insertEntity */
      static void insertEntity ( const typename BaseType::template Codim< 0 >::Entity & )
      {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::removeEntity */
      static void removeEntity ( const typename BaseType::template Codim< 0 >::Entity & )
      {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::backup */
      void backup () const {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::restore */
      void restore () {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::write */
      template< class T >
      void write ( OutStreamInterface< T > &stream ) const
      {}

      /** \copydoc Dune::Fem::ConsecutiveIndexSet::read */
      template< class T >
      void read ( InStreamInterface< T > &stream )
      {}

      /** \copydoc Dune::Fem::AdaptiveIndexSet::numberOfHoles */
      static constexpr int numberOfHoles ( ... ) noexcept
      {
        return 0;
      }

      /** \copydoc Dune::Fem::AdaptiveIndexSet::oldIndex */
      static constexpr int oldIndex ( ... ) noexcept
      {
        return 0;
      }

      /** \copydoc Dune::Fem::AdaptiveIndexSet::newIndex */
      static constexpr int newIndex ( ... ) noexcept
      {
        return 0;
      }

    private:
      const IndexSet &indexSet () const { return indexSet_; }

      const IndexSet &indexSet_;
    };



    namespace Capabilities
    {

      template< class IndexSet >
      struct isConsecutiveIndexSet< NonAdaptiveIndexSet< IndexSet > >
      {
        static const bool v = false;
      };

      template< class IndexSet >
      struct isAdaptiveIndexSet< NonAdaptiveIndexSet< IndexSet > >
      {
        static const bool v = false;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSETWRAPPER_HH
