#ifndef DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH
#define DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH

#include <cassert>

#include <type_traits>

#include <dune/common/visibility.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    class HasLocalFunction;

  } // namespace Fem


  namespace FemPy
  {

    // getGridFunctionWrapper
    // ----------------------

    DUNE_EXPORT inline pybind11::object getGridFunctionWrapper ()
    {
      // store module as static variable. No need for Singleton here since
      // it can be loaded multiple times.
      static pybind11::object o = pybind11::module::import( "dune.ufl" ).attr( "GridFunction" );
      return o;
    }

  } // namespace FemPy

} // namespace Dune


namespace pybind11
{

  namespace detail
  {

    // type_caster for dune-fem grid functions
    // ---------------------------------------

    template< class T >
    struct type_caster< T, std::enable_if_t< std::is_base_of< Dune::Fem::HasLocalFunction, T >::value > >
      : public type_caster_base< T >
    {
      typedef type_caster_base< T > Base;

      bool load ( handle src, bool convert )
      {
        if( isinstance( src, Dune::FemPy::getGridFunctionWrapper() ) )
          return Base::load( getattr( src, "__impl__" ), convert );
        else
          return Base::load( src, convert );
      }

/*
      template< class V >
      static handle cast ( V &&v, return_value_policy policy, handle parent )
      {
        pybind11::handle obj = Base::cast( std::forward< V >( v ), policy, parent );
        if( obj )
        {
          return Dune::FemPy::getGridFunctionWrapper()(obj);
          tuple args( 1 );
          args[ 0 ] = reinterpret_steal< object >( obj );
          assert( args.ptr() );
          return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
        }
        else
          return obj;
      }
*/
    };

  } // namespace detail

} // namespace pybind11

#endif // #ifndef DUNE_FEMPY_PYBIND11_GRIDFUNCTION_HH
