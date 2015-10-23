#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH

//- dune-fem includes
#include <dune/fem/gridpart/common/capabilities.hh>

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
        // either implement this or leaf it away !!!
//        static const bool v = hasGrid< typename CoordFunction::GridPartType >::v;
        static const bool v = false;
      };


      template< class CoordFunction >
      class hasSingleGeometryType< GeoGridPart< CoordFunction > >
      {
        typedef typename CoordFunction::GridPartType HostGridPartType;

      public:
        static const bool v = hasSingleGeometryType< HostGridPartType >::v;
        static const unsigned int topologyId = hasSingleGeometryType< HostGridPartType >::topologyId;
      };


      template< class CoordFunction >
      struct isCartesian< GeoGridPart< CoordFunction > >
      {
        static const bool v = false;
      };


      template< class CoordFunction, int codim  >
      struct hasEntity< GeoGridPart< CoordFunction >, codim >
      {
        // disable codim > 0 && < dim entities because of missing interface for subIndex method
        // once this is implemented we can simply use hasEntity for the HostGridPart.
        static const bool v = ( codim == 0 || codim == CoordFunction::GridPartType :: dimension ) ?
            hasEntity< typename CoordFunction::GridPartType, codim >::v : false ;
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
