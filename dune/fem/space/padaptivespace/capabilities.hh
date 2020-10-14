#ifndef DUNE_FEM_SPACE_PADAPTIVESPACE_CAPABILITIES_HH
#define DUNE_FEM_SPACE_PADAPTIVESPACE_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/padaptivespace/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      // Dune::Fem::Capabilities for PAdaptiveDGSpace
      // --------------------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< PAdaptiveDGSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };



      // Dune::Fem::Capabilities for PAdaptiveLagrangeSpace
      // --------------------------------------------------

      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasFixedPolynomialOrder< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct hasStaticPolynomialOrder< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
        static const int order = polOrder;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isContinuous< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isLocalized< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct isAdaptive< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct threadSafe< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int polOrder, class Storage >
      struct viewThreadSafe< PAdaptiveLagrangeSpace< FunctionSpace, GridPart, polOrder, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_PADAPTIVESPACE_CAPABILITIES_HH
