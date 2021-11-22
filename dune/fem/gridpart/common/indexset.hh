#ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH
#define DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH

#include <type_traits>
#include <utility>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------
    template< class Traits >
    class IndexSet;
    template< class Traits >
    class ConsecutiveIndexSet;
    template< class Traits >
    class AdaptiveIndexSet;


    namespace Capabilities
    {

      /** \brief specialize with \b true if index set implements the
       *         dune-fem interface for index sets
       *
       *  \note default value is \b true if index set is derived from
       *        IndexSet
       */
      template< class IndexSet >
      struct isDuneFemIndexSet
      {
        template< class Traits >
        static std::true_type __isDuneFemIndexSet ( const Dune::Fem::IndexSet< Traits > & );

        static std::false_type __isDuneFemIndexSet ( ... );

      public:
        static const bool v = decltype( __isDuneFemIndexSet( std::declval< IndexSet >() ) )::value;
      };


      // isConsecutiveIndexSet
      // ---------------------

      /** \brief specialize with \b true if index set implements the
       *         interface for consecutive index sets
       *
       *  \note default value is \b true if index set is derived from
       *        ConsecutiveIndexSet
       */
      template< class IndexSet >
      struct isConsecutiveIndexSet
      {
        template< class Traits >
        static std::true_type __isConsecutive ( const ConsecutiveIndexSet< Traits > & );

        static std::false_type __isConsecutive ( ... );

      public:
        static const bool v = decltype( __isConsecutive( std::declval< IndexSet >() ) )::value;
      };


      // isAdaptiveIndexSet
      // ------------------

      /** \brief specialize with \b true if index set implements the
       *         interface for adaptive index sets
       *
       *  \note default value is \b true if index set is derived from
       *        AdaptiveIndexSet
       */
      template< class IndexSet >
      class isAdaptiveIndexSet
      {
        template< class Traits >
        static std::true_type __isAdaptive ( const AdaptiveIndexSet< Traits > & );

        static std::false_type __isAdaptive ( ... );

      public:
        static const bool v = decltype( __isAdaptive( std::declval< IndexSet >() ) )::value;
      };



#ifndef DOXYGEN

      template< class IndexSet >
      struct isConsecutiveIndexSet< const IndexSet >
        : public isConsecutiveIndexSet< IndexSet >
      {};

      template< class IndexSet >
      struct isAdaptiveIndexSet< const IndexSet >
        : public isAdaptiveIndexSet< IndexSet >
      {};

#endif // #ifndef DOXYGEN

    } // namespace Capabilites



    // IndexSet
    // --------

    /** \brief interface documentation for (grid part) index sets
     *
     *  \note This interface is a direct copy of index sets as used in
     *  dune-grid.
     */
    template< class Traits >
    class IndexSet
    {
    public:
      /** \brief grid dimension */
      static const int dimension = Traits::dimension;

      template< int codim >
      struct Codim
      {
        /** \brief entity type */
        typedef typename Traits::template Codim< codim >::Entity Entity;
      };

      /** \brief index type */
      typedef typename Traits::IndexType IndexType;

      /** \brief geometry type range type */
      typedef typename Traits::Types Types;

    protected:
      IndexSet () = default;

    public:
      /** \brief return range of geometry types */
      Types types ( int codim ) const
      {
        return impl().types( codim );
      }

      /** \brief return \b true if entity has index */
      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        return impl().contains( entity );
      }

      /** \brief return number of entities of given type */
      IndexType size ( GeometryType type ) const
      {
        return impl().size( type );
      }

      /** \brief return number of entities of given codimension */
      IndexType size ( int codim ) const
      {
        return impl().size( codim );
      }

      /** \brief return index for given entity */
      template< class Entity >
      IndexType index ( const Entity &entity ) const
      {
        return index< Entity::codimension >( entity );
      }

      /** \brief return index for given entity */
      template< int codim >
      IndexType index ( const typename Codim< codim >::Entity &entity ) const
      {
        return impl().template index< codim >( entity );
      }

      /** \brief return index for given subentity */
      template< class Entity >
      IndexType subIndex ( const Entity &entity, int i, unsigned int cd ) const
      {
        return subIndex< Entity::codimension >( entity, i, cd );
      }

      /** \brief return index for given subentity */
      template< int codim >
      IndexType subIndex ( const typename Codim< codim >::Entity &entity, int i, unsigned int cd ) const
      {
        return impl().template subIndex< codim >( entity, i, cd );
      }

      /** \brief receive request for codimension support in case index set is adaptive */
      void requestCodimensions ( const std::vector< int >& codimensions ) const
      {
      }

    protected:
      const typename Traits::IndexSetType &impl () const
      {
        return static_cast< const typename Traits::IndexSetType & >( *this );
      }
    };

    // ConsecutiveIndexSet
    // -------------------

    /** \brief extended interface for consecutive index sets
     *
     *  \note IndexSets implementing this extended interface can be managed by
     *        the DofManager
     */
    template< class Traits >
    class ConsecutiveIndexSet
      : public IndexSet< Traits >
    {
      typedef IndexSet< Traits > BaseType;

    protected:
      using BaseType::impl;

      ConsecutiveIndexSet () = default;

    public:
      /** \name Adaptation
       *  \{
       */

      /** \brief please doc me */
      void resize () { impl().resize(); }

      /** \brief please doc me */
      bool compress () { return impl().compress(); }

      /** \brief please doc me */
      void insertEntity ( const typename BaseType::template Codim< 0 >::Entity &entity )
      {
        impl().insertEntity( entity );
      }

      /** \brief please doc me */
      void removeEntity ( const typename BaseType::template Codim< 0 >::Entity &entity )
      {
        impl().removeEntity( entity );
      }

      /** \brief please doc me */
      void backup () const { impl().backup(); }

      /** \brief please doc me */
      void restore () { impl().restore(); }

      /** \brief please doc me */
      template< class T >
      void write ( OutStreamInterface< T > &stream ) const
      {
        impl().write( stream );
      }

      /** \brief please doc me */
      template< class T >
      void read ( InStreamInterface< T > &stream )
      {
        impl().read( stream );
      }

      /** \} */

    protected:
      typename Traits::IndexSetType &impl ()
      {
        const typename Traits::IndexSetType &impl = BaseType::impl();
        return const_cast< typename Traits::IndexSetType & >( impl );
      }
    };



    // AdaptiveIndexSet
    // ----------------

    /** \brief extended interface for adaptive, consecutive index sets
     *
     *  \note IndexSets implementing this extended interface can be used with
     *        index set based adaptive Dof mappers
     */
    template< class Traits >
    class AdaptiveIndexSet
      : public ConsecutiveIndexSet< Traits >
    {
      typedef ConsecutiveIndexSet< Traits > BaseType;

    protected:
      using BaseType::impl;

      AdaptiveIndexSet () = default;

    public:
      /** \name Adaptation
       *  \{
       */

      /** \brief return number of holes for given type */
      int numberOfHoles ( GeometryType type ) const
      {
        return impl().numberOfHoles( type );
      }

      /** \brief return old index for given hole and type */
      int oldIndex ( int hole, GeometryType type ) const
      {
        return impl().oldIndex( hole, type );
      }

      /** \brief return new index for given hole and type */
      int newIndex ( int hole, GeometryType type ) const
      {
        return impl().newIndex( hole, type );
      }

      /** \} */

    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH
