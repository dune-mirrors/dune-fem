#ifndef DUNE_FEM_SPACE_RANNACHERTUREK_CAPABILITIES_HH
#define DUNE_FEM_SPACE_RANNACHERTUREK_CAPABILITIES_HH

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/space/common/capabilities.hh>
#include <dune/fem/space/rannacherturek/declaration.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Capabilities
    {

      template< class FunctionSpace, class GridPart, class Storage >
      struct hasFixedPolynomialOrder< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct hasStaticPolynomialOrder< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
        static const int order = 1;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct isContinuous< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct isLocalized< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct isAdaptive< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct threadSafe< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = false;
      };


      template< class FunctionSpace, class GridPart, class Storage >
      struct viewThreadSafe< RannacherTurekDiscreteFunctionSpace< FunctionSpace, GridPart, Storage > >
      {
        static const bool v = true;
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_RANNACHERTUREK_CAPABILITIES_HH
