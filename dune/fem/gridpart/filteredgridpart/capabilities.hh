#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH

#include <dune/geometry/type.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class, class, bool >
    class FilteredGridPart;



    namespace GridPartCapabilities
    {

      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
      struct hasGrid< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > >
      {
        static const bool v = hasGrid< HostGridPartImp >::v;
      };


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
      struct hasSingleGeometryType< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > >
      {
        static const bool v = hasSingleGeometryType< HostGridPartImp >::v;
        static const unsigned int topologyId = hasSingleGeometryType< HostGridPartImp >::topologyId;
      };


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
      struct isCartesian< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > >
      {
        static const bool v = isCartesian< HostGridPartImp >::v;
      };


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet, int codim >
      struct hasEntity< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet >, codim >
      {
        static const bool v = hasEntity< HostGridPartImp, codim >::v;
      };


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet, int codim >
      struct canCommunicate< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet >, codim >
      {
        static const bool v = canCommunicate< HostGridPartImp, codim >::v;
      };


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
      struct isConforming< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > >
      {
        static const bool v = isConforming< HostGridPartImp >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH
