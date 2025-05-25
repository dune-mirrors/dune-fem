#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH

#include <cassert>

#include <type_traits>
#include <utility>

#include <dune/grid/common/intersectioniterator.hh>

#include <dune/fem/gridpart/filteredgridpart/intersection.hh>

namespace Dune
{

  namespace Fem
  {

    // FilteredGridPartIntersectionIterator
    // ------------------------------------

    template< class GridPartFamily >
    class FilteredGridPartIntersectionIterator
    {
      typedef FilteredGridPartIntersectionIterator< GridPartFamily > ThisType;

      typedef typename std::remove_const_t< GridPartFamily >::Filter FilterType;
      typedef typename std::remove_const_t< GridPartFamily >::ExtraData ExtraData;
      typedef typename std::remove_const_t< GridPartFamily >::HostGridPart::IntersectionIteratorType HostIteratorType;

      typedef FilteredGridPartIntersection< GridPartFamily > IntersectionImpl;

    public:
      typedef Dune::Intersection< GridPartFamily, IntersectionImpl > Intersection;

      FilteredGridPartIntersectionIterator () = default;

      FilteredGridPartIntersectionIterator ( ExtraData data, HostIteratorType hostIterator )
        : data_( data ), hostIterator_( std::move( hostIterator ) )
      {}

      Intersection dereference () const { return Intersection( IntersectionImpl( data(), *hostIterator_ ) ); }

      bool equals ( const ThisType &other ) const { return (hostIterator() == other.hostIterator()); }

      void increment () { ++hostIterator_; }

      const FilterType &filter () const { return data()->filter(); }

      const HostIteratorType &hostIterator () const { return hostIterator_; }
      HostIteratorType &hostIterator () { return hostIterator_; }

    protected:
      ExtraData data() const { assert( data_ ); return data_; }

      ExtraData data_;
      HostIteratorType hostIterator_;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH
