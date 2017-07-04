#ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_CAPABILITIES_HH
#define DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/brezzidouglasmarini/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct hasFixedPolynomialOrder< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int ord, template< class > class Storage >
      struct hasStaticPolynomialOrder< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, ord, Storage > >
      {
        static const bool v = true;
        static const int order = ord;
      };


      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct isContinuous< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct isLocalized< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct isAdaptive< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct threadSafe< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int order, template< class > class Storage >
      struct viewThreadSafe< BDMDiscreteFunctionSpace< FunctionSpace, GridPart, order, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_BREZZIDOUGLASMARINI_CAPABILITIES_HH
