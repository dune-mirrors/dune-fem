#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-geometry includes
#include <dune/geometry/genericgeometry/topologytypes.hh>

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/misc/selectvalue.hh>


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
        static const unsigned int topologyId
          = SelectUnsignedValue< v, hasSingleGeometryType< HostGridPartImp >::topologyId, ~0 >::Value;
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


      template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
      struct isParallel< FilteredGridPart< HostGridPartImp, FilterImp, useFilteredIndexSet > >
      {
        static const bool v = isParallel< HostGridPartImp >::v;
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
