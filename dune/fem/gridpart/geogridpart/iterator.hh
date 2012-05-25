#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH

#include <dune/geometry/referenceelements.hh>

#include <dune/grid/common/entityiterator.hh>

#include <dune/fem/gridpart/geogridpart/entitypointer.hh>

namespace Dune
{

  namespace Fem
  {

    // GeoIteratorTraits
    // -----------------

    template< int codim, PartitionIteratorType pitype, class GridFamily >
    struct GeoIteratorTraits
    : public GeoEntityPointerTraits< codim, GridFamily >
    {
      typedef typename remove_const< GridFamily >::type::Traits::HostGridPartType HostGridPartType;

      typedef typename HostGridPartType::template Codim< codim >::template Partition< pitype >::IteratorType HostIteratorType;
    };



    // GeoIterator
    // -----------
    
    template< int codim, PartitionIteratorType pitype, class GridFamily >
    class GeoIterator
    : public GeoEntityPointer< GeoIteratorTraits< codim, pitype, GridFamily > >
    {
      typedef GeoEntityPointer< GeoIteratorTraits< codim, pitype, GridFamily > > Base;

    protected:
      typedef typename Base::HostIteratorType HostIteratorType;

      using Base::hostIterator_;
      using Base::releaseEntity;

    public:
      typedef typename Base::CoordFunctionType CoordFunctionType;

      GeoIterator ( const CoordFunctionType &coordFunction, const HostIteratorType &hostIterator )
      : Base( coordFunction, hostIterator )
      {}
      
      void increment ()
      {
        ++hostIterator_;
        releaseEntity();
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_ITERATOR_HH
