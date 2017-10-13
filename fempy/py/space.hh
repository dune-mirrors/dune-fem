#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/corepy/common/fmatrix.hh>
#include <dune/corepy/common/fvector.hh>

#include <dune/fempy/pybind11/pybind11.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // registerSpaceConstructor
      // -------------------------
      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls, std::false_type )
      {}
      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls, std::true_type )
      {
        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;
        cls.def( "__init__", [] ( Space &self, pybind11::object gridView ) {
            new( &self ) Space( gridPart< GridView >( gridView ) );
          }, pybind11::keep_alive< 1, 2 >() );
      }
      template< class Space, class... options >
      void registerSpaceConstructor ( pybind11::class_< Space, options... > cls )
      {
        typedef typename Space::GridPartType GridPart;
        registerSpaceConstructor( cls, std::is_constructible< Space, GridPart& >() );
      }

      // registerSpace
      // -------------

      template< class Space, class... options >
      void registerSpace ( pybind11::module module, pybind11::class_< Space, options... > cls )
      {
        typedef typename Space::FunctionSpaceType::RangeFieldType RangeFieldType;
        if( !std::is_same< RangeFieldType, double >::value )
        {
          Dune::CorePy::registerFieldVector<RangeFieldType>(module, std::make_integer_sequence<int, 10>());
          Dune::CorePy::registerFieldMatrix<RangeFieldType>(module, std::make_integer_sequence<int, 5>());
        }

        typedef typename Space::GridPartType GridPart;
        typedef typename GridPart::GridViewType GridView;

        cls.def_property_readonly( "dimRange", [] ( Space & ) -> int { return Space::dimRange; } );
        cls.def_property_readonly( "grid", [] ( Space &self ) -> GridView { return static_cast< GridView >( self.gridPart() ); } );
        cls.def_property_readonly( "order", [] ( Space &self ) -> int { return self.order(); } );
        cls.def_property_readonly( "size", [] ( Space &self ) -> int { return self.size(); } );

        registerSpaceConstructor( cls );
      }



      // getSpace
      // --------

      template< class T >
      pybind11::object getSpace ( const T &obj, pybind11::handle nurse = pybind11::handle() )
      {
        typedef std::decay_t< decltype( std::declval< const T & >().space() ) > Space;

        const Space &space = obj.space();
        auto pySpace = pybind11::reinterpret_borrow< pybind11::object >( pybind11::detail::get_object_handle( &space, pybind11::detail::get_type_info( typeid( Space ) ) ) );
#if 0
        // disable returning spaces not defined on Python side for now
        // Bug: pybind11 seems to require at least a move constructor to cast a reference to a Python object
        if( !pySpace )
        {
          pySpace = pybind11::object( pybind11::cast( space, pybind11::return_value_policy::reference ), true );
          if( !nurse )
            nurse = pybind11::detail::get_object_handle( &obj, pybind11::detail::get_type_info( typeid( T ) ) );
          pybind11::detail::keep_alive_impl( pySpace, nurse );
        }
#endif
        return pySpace;
      }

    } // namespace detail



    // registerSpace
    // -------------

    template< class Space, class... options >
    void registerSpace ( pybind11::module module, pybind11::class_< Space, options... > cls )
    {
      detail::registerSpace( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
