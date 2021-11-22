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

      template< class Entity >
      GeoIntersectionIterator ( const Entity &inside,
                                HostIntersectionIteratorType hostIterator )
      : coordFunction_( &inside.impl().coordFunction() ),
        insideGeo_( inside.geometry() ),
        hostIterator_( std::move( hostIterator ) )
      {}

      GeoIntersectionIterator ( const GeoIntersectionIterator& other )
      : coordFunction_( &other.coordFunction() ),
        insideGeo_( other.insideGeo_ ),
        hostIterator_( other.hostIterator_ )
      {}

      GeoIntersectionIterator ( )
      : coordFunction_( nullptr ),
        hostIterator_()
      {}

      GeoIntersectionIterator& operator= ( const GeoIntersectionIterator& other )
      {
        coordFunction_ = other.coordFunction_;
        insideGeo_.reset();
        insideGeo_.emplace( *other.insideGeo_ );
        hostIterator_ = other.hostIterator_ ;
        return *this;
      }

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
        assert( insideGeo_ );
        return IntersectionImplType( coordFunction(), *insideGeo_, *hostIterator_ );
      }

      const CoordFunctionType &coordFunction () const
      {
        assert( coordFunction_ );
        return *coordFunction_;
      }

    private:
      const CoordFunctionType *coordFunction_ = nullptr;
      std::optional< ElementGeometryType > insideGeo_;
      HostIntersectionIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_INTERSECTIONITERATOR_HH
