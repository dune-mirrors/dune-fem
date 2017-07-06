#ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH
#define DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class LFEMap, class FunctionSpace, template< class > class Storage = CachingStorage >
    class LocalFiniteElementSpace;



    namespace Capabilities
    {

      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct hasFixedPolynomialOrder< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct hasStaticPolynomialOrder< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
        static const int order = 111;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isContinuous< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isLocalized< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = true;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isAdaptive< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct threadSafe< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct viewThreadSafe< LocalFiniteElementSpace< LFEMap, FunctionSpace, Storage > >
      {
        static const bool v = false;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH
