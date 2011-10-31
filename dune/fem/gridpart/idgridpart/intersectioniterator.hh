#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/fem/gridpart/idgridpart/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIntersectionIterator
    // ----------------------

    template< class GridFamily >
    class IdIntersectionIterator
    {
      typedef IdIntersectionIterator< GridFamily > ThisType;

      typedef typename remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType::IntersectionIteratorType HostIntersectionIteratorType;

      typedef IdIntersection< const GridFamily > IntersectionImplType;

    public:
      typedef Dune::Intersection< const GridFamily, IdIntersection > Intersection;

      IdIntersectionIterator ( const HostIntersectionIteratorType &hostIterator )
      : intersection_( IntersectionImplType() ),
        hostIterator_( hostIterator )
      {}

      IdIntersectionIterator ( const ThisType &other )
      : intersection_( IntersectionImplType() ),
        hostIterator_( other.hostIterator_ )
      {}

      const ThisType &operator= ( const ThisType &other )
      {
        intersectionImpl() = IntersectionImplType();
        hostIterator_ = other.hostIterator_;
        return *this;
      }

      bool equals ( const ThisType &other ) const
      {
        return (hostIterator_ == other.hostIterator_);
      }
      
      void increment ()
      {
        ++hostIterator_;
        intersectionImpl() = IntersectionImplType();
      }

      const Intersection &dereference () const
      {
        if( !intersectionImpl() )
          intersectionImpl() = IntersectionImplType( *hostIterator_ );
        return intersection_;
      }

    private:
      IntersectionImplType &intersectionImpl () const
      {
        return intersection_.impl();
      }

      mutable Intersection intersection_;
      HostIntersectionIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH
