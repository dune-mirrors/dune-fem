#include <config.h>

#include <string>
#include <tuple>

#include <dune/common/fvector.hh>
#include <dune/common/std/utility.hh>

#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/operators.h>

template< int size >
static void registerFieldVector ( pybind11::module module, std::integral_constant< int, size > )
{
  typedef Dune::FieldVector< double, size > FV;

  static const std::string clsName = "FieldVector" + std::to_string( size );
  pybind11::class_< FV > cls( module, clsName.c_str() );

  cls.def( pybind11::init<>() );

  cls.def( "__getitem__", [] ( const FV &x, std::size_t i ) -> double {
      if( i < size )
        return x[ i ];
      else
        throw pybind11::index_error();
    } );

  cls.def( "__setitem__", [] ( FV &x, std::size_t i, double v ) {
      if( i < size )
        x[ i ] = v;
      else
        throw pybind11::index_error();
    } );

  cls.def( "__len__", [] ( const FV &x ) -> std::size_t { return size; } );

  cls.def( pybind11::self + pybind11::self );
  cls.def( pybind11::self += pybind11::self );
  cls.def( pybind11::self - pybind11::self );
  cls.def( pybind11::self -= pybind11::self );
  cls.def( pybind11::self *= double() );
  cls.def( pybind11::self /= double() );

  cls.def( "__repr__", [] ( const FV &x ) {
      std::string repr = "DUNE FieldVector: (";
      for( int i = 0; i < size; ++i )
        repr += (i > 0 ? ", " : "") + std::to_string( x[ i ] );
      repr += ")";
      return repr;
    } );

  cls.def_property_readonly( "one_norm", [] ( const FV &x ) { return x.one_norm(); } );
  cls.def_property_readonly( "two_norm", [] ( const FV &x ) { return x.two_norm(); } );
  cls.def_property_readonly( "two_norm2", [] ( const FV &x ) { return x.two_norm2(); } );
  cls.def_property_readonly( "infinity_norm", [] ( const FV &x ) { return x.infinity_norm(); } );
}

template< int... size >
static void registerFieldVector ( pybind11::module module, Dune::Std::integer_sequence< int, size... > )
{
  std::ignore = std::make_tuple( (registerFieldVector( module, std::integral_constant< int, size >() ), size)... );
}


PYBIND11_PLUGIN( common )
{
  pybind11::module module( "common" );

  registerFieldVector( module, Dune::Std::make_integer_sequence< int, 10 >() );

  return module.ptr();
}
