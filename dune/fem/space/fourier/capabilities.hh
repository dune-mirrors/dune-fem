#ifndef DUNE_FEM_SPACE_FOURIER_CAPABILITIES_HH
#define DUNE_FEM_SPACE_FOURIER_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/fourier/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, int Order >
      struct hasFixedPolynomialOrder< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct hasStaticPolynomialOrder< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = false;
        static const int order = -1;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct isContinuous< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct isLocalized< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct isAdaptive< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct threadSafe< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, int Order >
      struct viewThreadSafe< FourierDiscreteFunctionSpace< FunctionSpace, GridPart, Order > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_FOURIER_CAPABILITIES_HH
