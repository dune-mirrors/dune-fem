#ifndef DUNE_FEM_GRIDPART_CAPABILITIES_HH
#define DUNE_FEM_GRIDPART_CAPABILITIES_HH

//- dune-grid includes
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridenums.hh>


namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class > class GridPartDefault;
    template< class, PartitionIteratorType, bool > class AdaptiveLeafGridPart;
    template< class > class HierarchicGridPart;
    template< class > class LeafGridPart;
    template< class > class LevelGridPart;



    namespace GridPartCapabilities
    {
      /** \brief specialize with 'false' if grid part has no 
       *         underlying dune grid (default=true)
       */
      template< class GridPartType >
      struct hasGrid
      {
        static const bool v = true;
      };


      /** \brief specialize with 'true' for if the codimension 0 entity 
       *         of the grid part has only one possible geometry type
       *         (default=false, topologyid=undefined)
       */
      template< class GridPartType >
      struct hasSingleGeometryType
      {
        static const bool v = false;
        static const unsigned int topologyId = ~0u;
      };


      /** \brief specialize with 'true' if the grid part is 
       *         cartesian (default=false)
       */
      template< class GridPartType >
      struct isCartesian
      {
        static const bool v = false;
      };


      /** \brief specialize with 'true' for all codims that a 
                 grid implements entities for (default=false)
      */
      template< class GridPartType, int codim  >
      struct hasEntity
      {
        static const bool v = false;
      };

      
      /** \brief specialize with 'true' if implementation supports 
       *         parallelism (default=false)
       */
      template< class GridPartType >
      struct isParallel
      {
        static const bool v = false;
      };


      /** \brief specialize with 'true' for all codims that a 
       *         grid can communicate data on (default=false)
       */
      template< class GridPartType, int codim >
      struct canCommunicate
      {
        static const bool v = false;
      };


      /** \brief specialize with 'true' if implementation guarantees 
       *         conforming level grids. (default=false)
       */
      template< class GridPartType >
      struct isConforming
      {
        static const bool v = false;
      };



      /*
       * forward
       *   GridPartCapabilities::Something< const GridPartType >
       * to
       *   GridPartCapabilities::Something< GridPartType >
       */

      template< class GridPartType >
      struct hasGrid< const GridPartType >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::hasGrid< GridPartType >::v;
      };


      template< class GridPartType >
      struct hasSingleGeometryType< const GridPartType >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::v;
        static const unsigned int topologyId 
          = Dune::Fem::GridPartCapabilities::hasSingleGeometryType< GridPartType >::topologyId;;
      };


      template< class GridPartType >
      struct isCartesian< const GridPartType >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isCartesian< GridPartType >::v;
      };


      template< class GridPartType, int codim  >
      struct hasEntity< const GridPartType, codim >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::hasEntity< GridPartType, codim >::v;
      };


      template< class GridPartType >
      struct isParallel< const GridPartType >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPartType >::v;
      };


      template< class GridPartType, int codim >
      struct canCommunicate< const GridPartType, codim >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::canCommunicate< GridPartType, codim >::v;
      };


      template< class GridPartType >
      struct isConforming< const GridPartType >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isConforming< GridPartType >::v;
      };

    } // end namespace GridPartCapabilities

  } // end namespace Fem

} // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_CAPABILITIES_HH
