#ifndef DUNE_FEM_SPACE_LAGRANGE_CAPABILITIES_HH
#define DUNE_FEM_SPACE_LAGRANGE_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/lagrange/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {
      // LagrangeDiscreteFunctionSpace
      //------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isConforming< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };



      // DynamicLagrangeDiscreteFunctionSpace
      //-------------------------------------

      template< class FunctionSpace, class GridPart, class Storage >
      struct hasFixedPolynomialOrder< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      // hasStaticPolynomialOrder we use default

      template< class FunctionSpace, class GridPart, class Storage >
      struct isContinuous< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isConforming< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct isLocalized< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct isAdaptive< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct threadSafe< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct viewThreadSafe< DynamicLagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_CAPABILITIES_HH
