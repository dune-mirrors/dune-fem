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
      typedef Dune::Intersection< const GridFamily, IntersectionImplType > Intersection;
      typedef typename Traits::ExtraData  ExtraData;

      IdIntersectionIterator ( ExtraData data, const HostIntersectionIteratorType &hostIterator )
      : intersection_( IntersectionImplType( data ) ),
        hostIterator_( hostIterator )
      {}

      IdIntersectionIterator ( const ThisType &other )
      : intersection_( IntersectionImplType( other.data() ) ),
        hostIterator_( other.hostIterator_ )
      {}

      const ThisType &operator= ( const ThisType &other )
      {
        intersectionImpl() = IntersectionImplType( other.data() );
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
        intersectionImpl() = IntersectionImplType( data() );
      }

      const Intersection &dereference () const
      {
        if( !intersectionImpl() )
          intersectionImpl() = IntersectionImplType( data(), *hostIterator_ );
        return intersection_;
      }

    protected:
      IntersectionImplType &intersectionImpl () const
      {
        return intersection_.impl();
      }

      ExtraData data () const { return intersectionImpl().data(); }

    protected:
      mutable Intersection intersection_;
      HostIntersectionIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH
