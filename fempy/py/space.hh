#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/corepy/common/fmatrix.hh>
#include <dune/corepy/common/fvector.hh>

#include <dune/corepy/pybind11/complex.h>
#include <dune/corepy/pybind11/pybind11.h>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerSpace
      // -------------

      template< class Space, class Cls >
      void registerSpace ( pybind11::module module, Cls &cls )
      {
        typedef typename Space::FunctionSpaceType::RangeFieldType RangeFieldType;
        if( !std::is_same< RangeFieldType, double >::value )
        {
          Dune::CorePy::registerFieldVector<RangeFieldType>(module, std::make_integer_sequence<int, 10>());
          Dune::CorePy::registerFieldMatrix<RangeFieldType>(module, std::make_integer_sequence<int, 5>());
        }

        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;

        cls.def_property_readonly( "grid", [] ( Space &spc ) -> GridView { return static_cast< GridView >( spc.gridPart() ); } );
        cls.def_property_readonly( "order", [] ( Space &spc ) -> int { return spc.order(); } );
        cls.def_property_readonly( "size", [] ( Space &spc ) -> int { return spc.size(); } );

        cls.def( "__init__", [] ( Space &instance, pybind11::object gridView ) {
            new( &instance ) Space( gridPart< GridView >( gridView ) );
          }, pybind11::keep_alive< 1, 2 >() );
      }

    } // namespace detail



    // registerSpace
    // -------------

    template< class Space, class Holder, class AliasType >
    void registerSpace ( pybind11::module module, pybind11::class_<Space,Holder,AliasType> &cls )
    {
      detail::registerSpace<Space>(module,cls);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
