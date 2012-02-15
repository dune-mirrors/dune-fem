#ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>

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
        static const bool v = true; 
      };


      template< class CoordFunction >
      struct hasSingleGeometryType< GeoGridPart< CoordFunction > >
      {
        static const bool v 
          = Dune::Capabilities::hasSingleGeometryType< typename CoordFunction::GridType >::v;
        static const unsigned int topologyId
          = Dune::Capabilities::hasSingleGeometryType< typename CoordFunction::GridType >::topologyId;
      };


      template< class CoordFunction >
      struct isCartesian< GeoGridPart< CoordFunction > >
      {
        static const bool v = false; 
      };


      template< class CoordFunction, int codim  >
      struct hasEntity< GeoGridPart< CoordFunction >, codim >
      {
        static const bool v = hasEntity< typename CoordFunction::GridType, codim >::v; 
      };


      template< class CoordFunction >
      struct isParallel< GeoGridPart< CoordFunction > >
      {
        static const bool v = isParallel< typename CoordFunction::GridType >::v;
      };


      template< class CoordFunction, int codim >
      struct canCommunicate< GeoGridPart< CoordFunction >, codim >
      {
        static const bool v = canCommunicate< typename CoordFunction::GridType, codim >::v;
      };


      template< class CoordFunction >
      struct isConforming< GeoGridPart< CoordFunction > >
      {
        static const bool v = isConforming< typename CoordFunction::GridType >::v;
      };

    } // end namespace GridPartCapabilities

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_GEOGRIDPART_CAPABILITIES_HH
