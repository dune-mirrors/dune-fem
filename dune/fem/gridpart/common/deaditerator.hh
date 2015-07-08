#ifndef DUNE_FEM_GRIDPART_COMMON_DEADINTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_COMMON_DEADINTERSECTIONITERATOR_HH

#include <type_traits>

#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/intersection.hh>
#include <dune/grid/common/intersectioniterator.hh>

namespace Dune
{

  namespace Fem
  {

    // DeadIterator
    // ------------

    template< class E >
    struct DeadIterator
    {
      typedef E Entity;

      static const int codimension = Entity::codimension;

      void increment ()
      {
        DUNE_THROW( InvalidStateException, "Trying to increment a dead iterator." );
      }

      Entity dereference () const
      {
        DUNE_THROW( InvalidStateException, "Trying to dereference a dead iterator." );
      }

      bool equals ( const DeadIterator & ) const
      {
        DUNE_THROW( InvalidStateException, "Trying to compare a dead iterator." );
      }
    };



    // DeadIntersection
    // ----------------

    template< class GridFamily >
    class DeadIntersection
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

    public:
      typedef typename std::remove_const< GridFamily >::type::ctype ctype;

      static const int dimension = std::remove_const< GridFamily >::type::dimension;
      static const int dimensionworld = std::remove_const< GridFamily >::type::dimensionworld;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 1 >::Geometry Geometry;
      typedef typename Traits::template Codim< 1 >::LocalGeometry LocalGeometry;

      Entity inside () const
      {
        DUNE_THROW( InvalidStateException, "Call to inside on dead intersection." );
      }

      Entity outside () const
      {
        DUNE_THROW( InvalidStateException, "Call to outside on dead intersection." );
      }

      bool boundary () const
      {
        DUNE_THROW( InvalidStateException, "Call to boundary on dead intersection." );
      }

      bool conforming () const
      {
        DUNE_THROW( InvalidStateException, "Call to conforming on dead intersection." );
      }

      bool neighbor () const
      {
        DUNE_THROW( InvalidStateException, "Call to neighbor on dead intersection." );
      }

      int boundaryId () const
      {
        DUNE_THROW( InvalidStateException, "Call to boundaryId on dead intersection." );
      }

      size_t boundarySegmentIndex () const
      {
        DUNE_THROW( InvalidStateException, "Call to boundarySegmentIndex on dead intersection." );
      }

      const LocalGeometry &geometryInInside () const
      {
        DUNE_THROW( InvalidStateException, "Call to geometryInInside on dead intersection." );
      }

      const LocalGeometry &geometryInOutside () const
      {
        DUNE_THROW( InvalidStateException, "Call to geometryInOutside on dead intersection." );
      }

      const Geometry &geometry () const
      {
        DUNE_THROW( InvalidStateException, "Call to geometry on dead intersection." );
      }

      GeometryType type () const
      {
        DUNE_THROW( InvalidStateException, "Call to type on dead intersection." );
      }

      int indexInInside () const
      {
        DUNE_THROW( InvalidStateException, "Call to indexInInside on dead intersection." );
      }

      int indexInOutside () const
      {
        DUNE_THROW( InvalidStateException, "Call to indexInOutside on dead intersection." );
      }

      FieldVector< ctype, dimensionworld >
      integrationOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        DUNE_THROW( InvalidStateException, "Call to integrationOuterNormal on dead intersection." );
      }

      FieldVector< ctype, dimensionworld >
      outerNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        DUNE_THROW( InvalidStateException, "Call to outerNormal on dead intersection." );
      }

      FieldVector< ctype, dimensionworld >
      unitOuterNormal ( const FieldVector< ctype, dimension-1 > &local ) const
      {
        DUNE_THROW( InvalidStateException, "Call to unitOuterNormal on dead intersection." );
      }

      FieldVector< ctype, dimensionworld > centerUnitOuterNormal () const
      {
        DUNE_THROW( InvalidStateException, "Call to centerUnitOuterNormal on dead intersection." );
      }
    };



    // DeadIntersectionIterator
    // ------------------------

    template< class GridFamily >
    class DeadIntersectionIterator
    {
      typedef DeadIntersectionIterator< GridFamily > ThisType;

      typedef DeadIntersection< const GridFamily > DeadIntersectionType;

    public:
      typedef Dune::Intersection< const GridFamily, DeadIntersectionType > Intersection;

      DeadIntersectionIterator ()
      {}

      bool equals ( const ThisType &other ) const
      {
        return true;
      }

      void increment ()
      {
        DUNE_THROW( InvalidStateException, "Trying to increment a dead intersection iterator." );
      }

      const Intersection &dereference () const
      {
        DUNE_THROW( InvalidStateException, "Trying to dereference a dead intersection iterator." );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_DEADINTERSECTIONITERATOR_HH
