#ifndef DUNE_FEMPY_PYBIND11_SPACE_HH
#define DUNE_FEMPY_PYBIND11_SPACE_HH

#include <cassert>

#include <type_traits>

#include <dune/common/visibility.hh>

#include <dune/fem/space/common/discretefunctionspace.hh>

#include <dune/python/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // getGridSpaceWrapper
    // ----------------------

    DUNE_EXPORT inline pybind11::object getSpaceWrapper ()
    {
      static pybind11::object o = pybind11::module::import( "dune.ufl" ).attr( "FemSpace" );
      return o;
    }

  } // namespace FemPy

} // namespace Dune


namespace pybind11
{

  namespace detail
  {

    // type_caster for dune-fem discrete function spaces
    // -------------------------------------------------

    template< class T >
    struct type_caster< T, std::enable_if_t<
      std::is_same< Dune::Fem::DFSpaceIdentifier,
                    std::decay_t< decltype( std::declval<T>().type() ) > >::value > >
      : public type_caster_base< T >
    {
      typedef type_caster_base< T > Base;

      bool load ( handle src, bool convert )
      {
        if( isinstance( src, Dune::FemPy::getSpaceWrapper() ) )
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
          // return Dune::FemPy::getSpaceWrapper()(obj);
          tuple args( 1 );
          args[ 0 ] = reinterpret_steal< object >( obj );
          assert( args.ptr() );
          return PyObject_Call( Dune::FemPy::getSpaceWrapper().ptr(), args.ptr(), nullptr );
        }
        else
          return obj;
      }
      */
    };
  } // namespace detail

} // namespace pybind11
#endif // #ifndef DUNE_FEMPY_PYBIND11_SPACE_HH
