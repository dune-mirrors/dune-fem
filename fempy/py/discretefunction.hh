#ifndef DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#define DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/extensions.h>

#include <dune/fem/space/common/interpolate.hh>
#include <dune/fempy/py/function/grid.hh>
#include <dune/fempy/py/grid/function.hh>

namespace Dune
{

  namespace FemPy
  {

    // interpolant
    // -----------
    template< class DF, class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
    void interpolant ( DF &df, std::string name, const GF &gf )
    {
      Fem::interpolate( gf, df );
    }

    template< class DF >
    void interpolant ( DF &df, typename DF::RangeType value )
    {
      return interpolant( df, simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 ) );
    }


    // registerDiscreteFunction
    // -------------

    template< class DF >
    void registerDiscreteFunction ( pybind11::module module )
    {
      typedef typename DF::DiscreteFunctionSpaceType Space;
      typedef typename DF::GridPartType GridPart;
      typedef typename DF::RangeType Value;

      auto cls = detail::registerGridFunction< DF >( module, "DiscreteFunction" );

      detail::clsVirtualizedGridFunction< GridPart, Value >( module ).def( "__init__", [] ( VirtualizedGridFunction< GridPart, Value > &instance, DF &df ) {
          new (&instance) VirtualizedGridFunction< GridPart, Value >( pyGridFunction( df ) );
        } );
      pybind11::implicitly_convertible< DF, VirtualizedGridFunction< GridPart, Value > >();

      cls.def_property_readonly( "space", [](DF &df) -> const typename DF::DiscreteFunctionSpaceType& {return df.space();} );
      cls.def("clear", [] (DF &instance) { instance.clear(); } );

      cls.def( "__init__", [] ( DF &instance, Space &space, const std::string &name ) {
          new( &instance ) DF( std::move(name), space );
        }, pybind11::keep_alive< 1, 2 >() );

      cls.def( "assign", [] ( DF &instance, const DF &other ) { instance.assign(other); } );

      typedef VirtualizedGridFunction< GridPart, typename Space::RangeType > GridFunction;
      cls.def( "interpolate", [] ( DF &df, const GridFunction &gf ) {
          Fem::interpolate( gf, df );
        } );
      cls.def( "interpolate", [] ( DF &df, typename Space::RangeType value ) {
          Fem::interpolate(
                simpleGridFunction( df.space().gridPart(), [ value ] ( typename DF::DomainType ) { return value; }
          ), df );
        } );


      typedef typename DF::DofVectorType DofVector;
      if (!pybind11::already_registered<DofVector>())
      {
        auto clsDof = pybind11::class_<DofVector>( module, "DofVector");
      }
      cls.def( "dofVector", [] ( DF &instance ) { return instance.dofVector(); } );
      cls.def( "assign", [] ( DF &instance, const DofVector &other ) { instance.dofVector() = other; } );

    }

  } // namespace FemPy

} // namespace Dune


#endif // DUNE_FEMPY_PY_DISCRETEFUNCTION_HH
