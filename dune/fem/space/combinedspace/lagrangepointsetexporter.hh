#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_LAGRANGEPOINTSETEXPORTER_HH
#define DUNE_FEM_SPACE_COMBINEDSPACE_LAGRANGEPOINTSETEXPORTER_HH


namespace Dune
{

  namespace Fem
  {

    namespace CombinedSpaceHelper
    {

      // forward declaration of LagrangeDiscreteFunctionSpace
      // ----------------------------------------------------
      template< class FunctionSpaceImp,
                class GridPartImp,
                int polOrder,
                class BasisFunctionStorageImp >
      class LagrangeDiscreteFunctionSpace;



      // Helper class for combined space to extract the Lagrange Point Set
      // for one of the contained spaces.
      // -----------------------------------------------------------------

      template <class DFSpace>
      struct LagrangePointSetExporter
      {
        // for arbitrary spaces nothing needs to be done
        LagrangePointSetExporter( const DFSpace& spc ) {}
      };


      // specialization for Lagrange spaces
      // ----------------------------------

      template <class FunctionSpaceImp, class GridPartImp, int polOrder, class BaseFunctionStorageImp >
      struct LagrangePointSetExporter<
          LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, BaseFunctionStorageImp > >
      {
        // export the tpye of LagrangePointSet
        typedef LagrangeDiscreteFunctionSpace< FunctionSpaceImp, GridPartImp, polOrder, BaseFunctionStorageImp > LagrangeSpaceType;
        typedef typename LagrangeSpaceType :: LagrangePointSetType LagrangePointSetType;

        LagrangePointSetExporter( const LagrangeSpaceType& spc )
        : lagrangeSpace_( spc ) {}

        // return the lagrangepointset for a given entity
        template <class Entity>
        const LagrangePointSetType& lagrangePointSet( const Entity& entity ) const
        {
          return lagrangeSpace_.lagrangePointSet( entity );
        }

       private:
        const LagrangeSpaceType& lagrangeSpace_;
      };

    } // namespace CombinedSpaceHelper

  } // namespace Fem

} // namespace Dune

#endif //#ifndef DUNE_FEM_SPACE_COMBINEDSPACE_LAGRANGEPOINTSETEXPORTER_HH
