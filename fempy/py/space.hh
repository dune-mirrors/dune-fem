#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/corepy/common/fmatrix.hh>
#include <dune/corepy/common/fvector.hh>

#include <dune/corepy/pybind11/complex.h>
#include <dune/corepy/pybind11/pybind11.h>

#include <dune/fempy/function/virtualizedgridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // registerSpace
    // -------------

    namespace detail
    {
      template< class Space, class Cls >
      void registerSpace ( pybind11::module module, Cls &cls )
      {
        typedef typename Space::GridPartType GridPartType;
        typedef typename Space::FunctionSpaceType::RangeFieldType RangeFieldType;
        if (!std::is_same<RangeFieldType,double>::value)
        {
          Dune::CorePy::registerFieldVector<RangeFieldType>(module, std::make_integer_sequence<int, 10>());
          Dune::CorePy::registerFieldMatrix<RangeFieldType>(module, std::make_integer_sequence<int, 5>());
        }

        typedef typename Space::GridPartType GridPart;

        cls.def_property_readonly( "grid", [](Space &sp) -> const GridPart& {return sp.gridPart();} );

        cls.def( "__init__", [] ( Space &instance, GridPart &grid ) {
            new( &instance ) Space( grid );
          }, pybind11::keep_alive< 1, 2 >() );
      }
    }
    template< class Space, class Holder, class AliasType >
    void registerSpace ( pybind11::module module, pybind11::class_<Space,Holder,AliasType> &cls )
    {
      detail::registerSpace<Space>(module,cls);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
