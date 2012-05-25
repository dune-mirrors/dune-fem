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

    template< class >
    class IdGridPart;



    namespace GridPartCapabilities
    {

      template< class HostGridPartType >
      struct hasGrid< IdGridPart< HostGridPartType > >
      {
        static const bool v = true; 
      };


      template< class HostGridPartType >
      struct hasSingleGeometryType< IdGridPart< HostGridPartType > >
      {
        static const bool v = hasSingleGeometryType< HostGridPartType >::v;
        static const unsigned int topologyId 
          = SelectUnsignedValue< v, hasSingleGeometryType< HostGridPartType >::topologyId, ~0 >::Value;
      };


      template< class HostGridPartType >
      struct isCartesian< IdGridPart< HostGridPartType > >
      {
        static const bool v = isCartesian< HostGridPartType >::v;
      };


      template< class HostGridPartType, int codim  >
      struct hasEntity< IdGridPart< HostGridPartType >, codim >
      {
        static const bool v = hasEntity< HostGridPartType, codim >::v; 
      };


      template< class HostGridPartType >
      struct isParallel< IdGridPart< HostGridPartType > >
      {
        static const bool v = isParallel< HostGridPartType >::v;
      };


      template< class HostGridPartType, int codim >
      struct canCommunicate< IdGridPart< HostGridPartType >, codim >
      {
        static const bool v = canCommunicate< HostGridPartType, codim >::v;
      };


      template< class HostGridPartType >
      struct isConforming< IdGridPart< HostGridPartType > >
      {
        static const bool v = isConforming< HostGridPartType >::v;
      };

    } // end namespace GridPartCapabilities

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_CAPABILITIES_HH
