#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/fem/gridpart/geometrygridpart/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartIntersectionIterator
    // ----------------------

    template< class GridFamily >
    class GeometryGridPartIntersectionIterator
    {
      typedef GeometryGridPartIntersectionIterator< GridFamily > ThisType;

      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType::IntersectionIteratorType HostIntersectionIteratorType;

      typedef GeometryGridPartIntersection< const GridFamily > IntersectionImplType;

    public:
      typedef Dune::Intersection< const GridFamily, IntersectionImplType > Intersection;

      template <class Entity>
      GeometryGridPartIntersectionIterator ( const Entity &inside,
             const HostIntersectionIteratorType &hostIterator )
      : hostIterator_( hostIterator ),
        intersection_( IntersectionImplType( inside.impl().gridFunction(), inside.geometry() ) )
      {}

      GeometryGridPartIntersectionIterator ( const ThisType &other )
      : hostIterator_( other.hostIterator_ ),
        intersection_( IntersectionImplType( other.intersection_.impl() ) )
      {}

      const ThisType &operator= ( const ThisType &other )
      {
        hostIterator_ = other.hostIterator_;
        intersection_.impl() = other.intersection_.impl();
        return *this;
      }

      bool equals ( const ThisType &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }

      void increment ()
      {
        ++hostIterator_;
        intersection_.impl().invalidate();
      }

      const Intersection &dereference () const
      {
        if( !intersection_.impl() )
          intersection_.impl().initialize( *hostIterator_ );
        return intersection_;
      }

    private:
      HostIntersectionIteratorType hostIterator_;
      mutable Intersection intersection_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_INTERSECTIONITERATOR_HH
