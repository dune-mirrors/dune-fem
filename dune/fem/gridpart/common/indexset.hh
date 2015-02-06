#ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH
#define DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH

#include <type_traits>
#include <vector>

#include <dune/common/deprecated.hh>

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
    class AdaptiveIndexSet;



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

      /** \brief return vector of geometry types used of given codimension */
      const std::vector< GeometryType > &geomTypes ( int codim ) const DUNE_DEPRECATED
      {
        return impl().geomTypes( codim );
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
      IndexType index ( const typename Codim< codim >::EntityType &entity ) const
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
      IndexType subIndex ( const typename Codim< codim >::EntityType &entity, int i, unsigned int cd ) const
      {
        return impl().template subIndex< codim >( entity, i, cd );
      }

    protected:
      const typename Traits::IndexSetType impl () const
      {
        return static_cast< const typename Traits::IndexSetType & >( *this );
      }
    };



    // AdaptiveIndexSet
    // ----------------

    /** \brief extended interface for adaptive, consecutive index sets
     */
    template< class Traits >
    class AdaptiveIndexSet
      : public IndexSet< Traits >
    {
      typedef IndexSet< Traits > BaseType;

    protected:
      using BaseType::impl;

      AdaptiveIndexSet () = default;

    public:
      /** \name Adaptation
       *  \{
       */

      /** \brief please doc me */
      bool consecutive () const { return impl().consecutive(); }

      /** \brief please doc me */
      void resize () { impl().resize(); }

      /** \brief please doc me */
      bool compress () { return impl().compress(); }

      /** \brief please doc me */
      template< class Entity >
      void insertEntity ( const Entity &entity )
      {
        insertEntity< Entity::codimension >( entity );
      }

      /** \brief please doc me */
      template< int codim >
      void insertEntity ( const typename BaseType::template Codim< codim >::EntityType &entity )
      {
        impl().template insertEntity< codim >( entity );
      }

      /** \brief please doc me */
      template< class Entity >
      void removeEntity ( const Entity &entity )
      {
        removeEntity< Entity::codimension >( entity );
      }

      /** \brief please doc me */
      template< int codim >
      void removeEntity ( const typename BaseType::template Codim< codim >::EntityType &entity )
      {
        impl().template removeEntity< codim >( entity );
      }

      /** \} */

      /** \name Persistency
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

      /** \name Input/Output
       *  \{
       */

      /** \brief please doc me */
      void backup () const { impl().backup(); }

      /** \brief please doc me */
      void restore () { impl().restore(); }

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

      /** \} */
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_INDEXSET_HH
