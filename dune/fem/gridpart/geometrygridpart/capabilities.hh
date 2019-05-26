#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_CAPABILITIES_HH

#include <dune/common/version.hh>

#include <dune/fem/gridpart/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class >
    class GeometryGridPart;



    namespace GridPartCapabilities
    {

      template< class GridFunctionType >
      struct hasGrid< GeometryGridPart< GridFunctionType > >
      {
        static const bool v = true;
      };


      template< class GridFunctionType >
      struct hasSingleGeometryType< GeometryGridPart< GridFunctionType > >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = hasSingleGeometryType< HostGridPartType >::v;
        static const unsigned int topologyId = hasSingleGeometryType< HostGridPartType >::topologyId;
      };


      template< class GridFunctionType >
      struct isCartesian< GeometryGridPart< GridFunctionType > >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = isCartesian< HostGridPartType >::v;
      };

/*
      template< class GridFunctionType  >
      struct hasEntity< GeometryGridPart< GridFunctionType >, 0 >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = hasEntity< HostGridPartType, 0>::v;
      };
*/
      template< class GridFunctionType, int codim  >
      struct hasEntity< GeometryGridPart< GridFunctionType >, codim >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = false;
      };


      template< class GridFunctionType, int codim  >
      struct canCommunicate< GeometryGridPart< GridFunctionType >, codim >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = canCommunicate< HostGridPartType, codim >::v;
      };


      template< class GridFunctionType >
      struct isConforming< GeometryGridPart< GridFunctionType > >
      {
        typedef typename GridFunctionType::GridPartType HostGridPartType;
        static const bool v = isConforming< HostGridPartType >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_CAPABILITIES_HH
