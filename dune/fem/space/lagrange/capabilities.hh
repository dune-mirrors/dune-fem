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

      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasFixedPolynomialOrder< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct hasStaticPolynomialOrder< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isContinuous< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = Dune::Fem::GridPartCapabilities::isConforming< GridPart >::v;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isLocalized< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct isAdaptive< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct threadSafe< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, template< class > class Storage >
      struct viewThreadSafe< LagrangeDiscreteFunctionSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_CAPABILITIES_HH
