#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH

#include <type_traits>

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/fem/gridpart/geometrygridpart/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartIntersectionIterator
    // ------------------------------------

    template< class GridFamily >
    class GeometryGridPartIntersectionIterator
    {
      typedef GeometryGridPartIntersectionIterator< GridFamily > ThisType;

      typedef typename std::remove_const_t< GridFamily >::Traits Traits;

      typedef typename Traits::HostGridPartType::IntersectionIteratorType HostIntersectionIteratorType;

      typedef typename Traits::template Codim< 0 >::Entity Entity;
      typedef typename Traits::template Codim< 0 >::Geometry ElementGeometry;

      typedef typename Traits::GridFunctionType GridFunctionType;

      typedef GeometryGridPartIntersection< const GridFamily > IntersectionImplType;

    public:
      typedef Dune::Intersection< const GridFamily, IntersectionImplType > Intersection;

      GeometryGridPartIntersectionIterator () = default;

      GeometryGridPartIntersectionIterator ( const Entity &inside, const HostIntersectionIteratorType &hostIterator )
        : hostIterator_( hostIterator ), gridFunction_( &inside.impl().gridFunction() ), insideGeo_( inside.geometry().impl() )
      {}

      bool equals ( const ThisType &other ) const { return (hostIterator_ == other.hostIterator_); }

      void increment () { ++hostIterator_; }

      Intersection dereference () const { return IntersectionImplType( *gridFunction_, insideGeo_, *hostIterator_ ); }

    private:
      HostIntersectionIteratorType hostIterator_;
      const GridFunctionType *gridFunction_ = nullptr;
      typename ElementGeometry::Implementation insideGeo_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH
