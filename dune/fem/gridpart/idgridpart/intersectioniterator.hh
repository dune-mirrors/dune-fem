#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH

#include <type_traits>
#include <utility>

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

      typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      typedef typename Traits::HostGridPartType::IntersectionIteratorType HostIntersectionIteratorType;

      typedef IdIntersection< const GridFamily > IntersectionImplType;

    public:
      typedef Dune::Intersection< const GridFamily, IntersectionImplType > Intersection;
      typedef typename Traits::ExtraData  ExtraData;

      IdIntersectionIterator () = default;

      IdIntersectionIterator ( ExtraData data, HostIntersectionIteratorType hostIterator )
      : data_( std::move( data ) ),
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
        return IntersectionImplType( data(), *hostIterator_ );
      }

      const ExtraData &data () const { return data_; }

    protected:
      ExtraData data_;
      HostIntersectionIteratorType hostIterator_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INTERSECTIONITERATOR_HH
