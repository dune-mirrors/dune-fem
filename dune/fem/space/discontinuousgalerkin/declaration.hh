#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DECLARATION_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DECLARATION_HH


namespace Dune
{

  namespace Fem
  {

    // DiscontinuousGalerkinSpace
    // --------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class DiscontinuousGalerkinSpace;


    // LagrangeDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LagrangeDiscontinuousGalerkinSpace;


    // LegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class LegendreDiscontinuousGalerkinSpace;


    // HierarchicLegendreDiscontinuousGalerkinSpace
    // ----------------------------------

    template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
    class HierarchicLegendreDiscontinuousGalerkinSpace;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_DECLARATION_HH
