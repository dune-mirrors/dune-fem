#ifndef DUNE_FEMPY_PY_REFERENCEELEMENTS_HH
#define DUNE_FEMPY_PY_REFERENCEELEMENTS_HH

#include <array>

#include <dune/geometry/referenceelements.hh>

#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerReferenceElement
    // ------------------------

    template< class ctype, int dim >
    auto registerReferenceElement_ ( pybind11::handle scope )
    {
      typedef typename Dune::ReferenceElement< ctype, dim > RefElement;

      static const std::string name = "ReferenceElement" + std::to_string( dim );
      pybind11::class_< RefElement > cls( scope, name.c_str());

      cls.def( "size", [ ]( const RefElement &e, int c ) { return e.size( c ); } );
      cls.def( "size", [ ]( const RefElement &e, int i, int c, int cc ) { return e.size( i, c, cc ); } );
      cls.def( "type", [ ]( const RefElement &e, int i, int c ) { return e.type( i, c ); } );

      cls.def( "subEntity", &RefElement::subEntity );
      cls.def( "position", &RefElement::position );
      cls.def_property_readonly( "center", [ ]( const RefElement &ref ) { return ref.position( 0, 0 ); } );
      cls.def( "volume", &RefElement::volume );
      cls.def( "integrationOuterNormal", &RefElement::integrationOuterNormal );

      cls.def_property_readonly( "type", [ ]( const RefElement &e ) { return e.type(); } );

      // registerGeometry<RefElement>(cls, std::make_integer_sequence<int, dim+1>());

      return cls;
    }

    template< class ctype, int... dim >
    void registerReferenceElement ( pybind11::handle scope, std::integer_sequence< int, dim ... >)
    {
      std::ignore = std::make_tuple( registerReferenceElement_< ctype, dim >( scope )... );
    }



    // registerReferenceElements
    // -------------------------

    template< class ctype, int dim >
    pybind11::object registerReferenceElements_ ( pybind11::handle scope )
    {
      typedef typename Dune::ReferenceElements< ctype, dim > RefElements;

      static const std::string name = "ReferenceElements" + std::to_string( dim );
      pybind11::class_< RefElements > cls( scope, name.c_str() );

      cls.def_static( "general", &RefElements::general, pybind11::return_value_policy::reference_internal );
      cls.def_static( "simplex", &RefElements::simplex, pybind11::return_value_policy::reference_internal );
      cls.def_static( "cube", &RefElements::cube, pybind11::return_value_policy::reference_internal );

      return pybind11::cast( RefElements() );
    }

    template< class ctype, int... dim >
    void registerReferenceElements ( pybind11::module module, std::integer_sequence< int, dim... > )
    {
      static const std::array< pybind11::object, sizeof...( dim ) > referenceElements = {{ registerReferenceElements_< ctype, dim >( module )... }};

      module.def( "ReferenceElements", []( int dimension ) { return referenceElements[ dimension ]; } );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_REFERENCEELEMENTS_HH
