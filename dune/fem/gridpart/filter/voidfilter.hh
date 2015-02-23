#ifndef DUNE_FEM_GRIDPART_FILTER_VOIDFILTER_HH
#define DUNE_FEM_GRIDPART_FILTER_VOIDFILTER_HH

namespace Dune
{

  namespace Fem
  {

    // VoidFilter
    // ------------

    template< class GridPart >
    class VoidFilter
    {
    public:
      typedef VoidFilter FilterType;

      typedef GridPart GridPartType;

      template < int cd >
      struct Codim
      {
        typedef typename GridPartType::template Codim< cd >::EntityType EntityType;
      };

      typedef typename Codim< 0 >::EntityType EntityType;

      template< class Entity >
      static bool contains ( const Entity & entity )
      {
        return true;
      }

      template< int cd >
      static bool contains ( const typename Codim< cd >::EntityType & entity )
      {
        return true;
      }

      template < class IntersectionIteratorType >
      static bool intersectionBoundary( const IntersectionIteratorType & it )
      {
        return true;
      }

      template < class IntersectionIteratorType >
      static int intersectionBoundaryId(const IntersectionIteratorType & it)
      {
        return 1;
      }

      template <class IntersectionIteratorType>
      static bool intersectionNeighbor( const IntersectionIteratorType & it )
      {
        return true;
      }

      template< class Intersection >
      static bool interiorIntersection( const Intersection & intersection )
      {
        return true;
      }

    }; // end RadialFilter

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTER_VOIDFILTER_HH
