#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTIONITERATOR_HH

#include <type_traits>
#include <utility>

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/fem/gridpart/geogridpart/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoIntersectionIterator
    // -----------------------

    template< class GridFamily >
    class GeoIntersectionIterator
    {
      typedef GeoIntersectionIterator< GridFamily > ThisType;

      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::CoordFunctionType CoordFunctionType;
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometryType;
      typedef typename Traits::HostGridPartType::IntersectionIteratorType HostIntersectionIteratorType;

      typedef GeoIntersection< const GridFamily > IntersectionImplType;

    public:
      typedef Dune::Intersection< const GridFamily, IntersectionImplType > Intersection;

      GeoIntersectionIterator () = default;

      template< class Entity >
      GeoIntersectionIterator ( const Entity &inside,
                                HostIntersectionIteratorType hostIterator )
      : coordFunction_( &inside.impl().coordFunction() ),
        insideGeo_( inside.geometry() ),
        hostIterator_( std::move( hostIterator ) )
      {}

      bool equals ( const ThisType &other ) const
      {
        return hostIterator_ == other.hostIterator_;
      }

      void increment ()
      {
        ++hostIterator_;
      }

      Intersection dereference () const
      {
        return IntersectionImplType( coordFunction(), insideGeo_, *hostIterator_ );
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      ElementGeometryType insideGeo_;
      HostIntersectionIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTIONITERATOR_HH
