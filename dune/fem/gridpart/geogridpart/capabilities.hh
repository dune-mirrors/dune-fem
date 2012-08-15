#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH

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
    class GeoGridPart;



    namespace GridPartCapabilities
    {

      template< class CoordFunction >
      struct hasGrid< GeoGridPart< CoordFunction > >
      {
        static const bool v = hasGrid< typename CoordFunction::GridPartType >::v;
      };


      template< class CoordFunction >
      class hasSingleGeometryType< GeoGridPart< CoordFunction > >
      {
        typedef typename CoordFunction::GridPartType HostGridPartType;

      public:
        static const bool v = hasSingleGeometryType< HostGridPartType >::v;
        static const unsigned int topologyId
          = SelectUnsignedValue< v, hasSingleGeometryType< HostGridPartType >::topologyId, ~0 >::Value;
      };


      template< class CoordFunction >
      struct isCartesian< GeoGridPart< CoordFunction > >
      {
        static const bool v = false; 
      };


      template< class CoordFunction, int codim  >
      struct hasEntity< GeoGridPart< CoordFunction >, codim >
      {
        static const bool v = hasEntity< typename CoordFunction::GridPartType, codim >::v; 
      };


      template< class CoordFunction >
      struct isParallel< GeoGridPart< CoordFunction > >
      {
        static const bool v = isParallel< typename CoordFunction::GridPartType >::v;
      };


      template< class CoordFunction, int codim >
      struct canCommunicate< GeoGridPart< CoordFunction >, codim >
      {
        static const bool v = canCommunicate< typename CoordFunction::GridPartType, codim >::v;
      };


      template< class CoordFunction >
      struct isConforming< GeoGridPart< CoordFunction > >
      {
        static const bool v = isConforming< typename CoordFunction::GridPartType >::v;
      };

    } // namespace GridPartCapabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH
