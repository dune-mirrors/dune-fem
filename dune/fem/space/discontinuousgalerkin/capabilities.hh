#ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_CAPABILITIES_HH
#define DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/discontinuousgalerkin/declaration.hh>

namespace Dune
{
 
  namespace Fem
  {

    namespace Capabilities
    {

      // Dune::Fem::Capabilities for DiscontinuousGalerkinSpace
      // ------------------------------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isParallel< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< DiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };



      // Dune::Fem::Capabilities for LagrangeDiscontinuousGalerkinSpace
      // --------------------------------------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isParallel< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< LagrangeDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };



      // Dune::Fem::Capabilities for LegendreDiscontinuousGalerkinSpace
      // --------------------------------------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isParallel< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isParallel< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< LegendreDiscontinuousGalerkinSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_DISCONTINUOUSGALERKIN_CAPABILITIES_HH
