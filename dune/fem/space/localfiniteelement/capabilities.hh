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

    template< class LFEMap, class FunctionSpace, template< class > class Storage >
    struct LocalFiniteElementDiscreteFunctionSpace;



    namespace Capabilities
    {

      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct hasFixedPolynomialOrder< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct hasStaticPolynomialOrder< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
        static const int order = 111;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isContinuous< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isLocalized< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = true;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct isAdaptive< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct threadSafe< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
      };


      template< class LFEMap, class FunctionSpace, template< class > class Storage >
      struct viewThreadSafe< LocalFiniteElementDiscreteFunctionSpace< LFEMap, FunctionSpace,  Storage > >
      {
        static const bool v = false;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LOCALFINITEELEMENT_CAPABILITIES_HH
