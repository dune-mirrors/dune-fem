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
      typedef typename std::remove_const_t< GridPartFamily >::HostGridPart::IntersectionIteratorType HostIteratorType;

      typedef FilteredGridPartIntersection< FilterType, typename HostIteratorType::Intersection > IntersectionImpl;

    public:
      typedef Dune::Intersection< GridPartFamily, IntersectionImpl > Intersection;

      FilteredGridPartIntersectionIterator () = default;

      FilteredGridPartIntersectionIterator ( const FilterType &filter, HostIteratorType hostIterator )
        : filter_( &filter ), hostIterator_( std::move( hostIterator ) )
      {}

      Intersection dereference () const { return Intersection( IntersectionImpl( filter(), *hostIterator_ ) ); }

      bool equals ( const ThisType &other ) const { return (hostIterator() == other.hostIterator()); }

      void increment () { ++hostIterator_; }

      const FilterType &filter () const { assert( filter_ ); return *filter_; }

      const HostIteratorType &hostIterator () const { return hostIterator_; }
      HostIteratorType &hostIterator () { return hostIterator_; }

    private:
      const FilterType *filter_ = nullptr;
      HostIteratorType hostIterator_;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_INTERSECTIONITERATOR_HH
