#ifndef DUNE_FEMPY_PY_SPACE_HH
#define DUNE_FEMPY_PY_SPACE_HH

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/function/vectorfunction.hh>
#include <dune/fem/space/common/interpolate.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/function/discrete.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // registerSpace
    // -------------

    template< class Space >
    void registerSpace ( pybind11::module module )
    {
      typedef typename Space::GridPartType GridPart;

      pybind11::class_< Space > cls( module, "Space" );

      cls.def( "__init__", [] ( Space &instance, GridPart &grid ) {
          new( &instance ) Space( grid );
        }, pybind11::keep_alive< 1, 2 >() );

      //typedef Fem::ManagedDiscreteFunction< Fem::VectorDiscreteFunction< Space, NumPyVector< double > > > DiscreteFunction;
      //typedef Fem::ManagedDiscreteFunction< Fem::VectorDiscreteFunction< Space, DynamicVector< double > > > DiscreteFunction;
      typedef Fem::AdaptiveDiscreteFunction< Space > DiscreteFunction;

      registerDiscreteFunction< DiscreteFunction >( module );

      typedef VirtualizedGridFunction< GridPart, typename Space::RangeType > GridFunction;
      cls.def( "interpolate", [] ( const Space &space, const GridFunction &gf ) {
          using Fem::interpolate;
          DiscreteFunction df( gf.name(), space );
          interpolate( gf, df );
          return df;
        }, pybind11::keep_alive< 0, 1 >() );

      cls.def( "interpolate", [] ( const Space &space, typename Space::RangeType value ) {
          using Fem::interpolate;
          DiscreteFunction df( "constant", space );
          interpolate( simpleGridFunction( space.gridPart(), [ value ] ( typename Space::DomainType ) { return value; }, 0 ), df );
          return df;
        }, pybind11::keep_alive< 0, 1 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
