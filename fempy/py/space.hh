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

    // interpolant
    // -----------

    template< class DF, class GF, std::enable_if_t< std::is_base_of< Fem::HasLocalFunction, GF >::value, int > = 0 >
    std::unique_ptr<DF> interpolant ( const typename DF::DiscreteFunctionSpaceType &space, std::string name, const GF &gf )
    {
      using Fem::interpolate;
      std::unique_ptr<DF> df( new DF(std::move( name ), space) );
      interpolate( gf, *df );
      return df;
    }

    template< class DF >
    std::unique_ptr<DF> interpolant ( const typename DF::DiscreteFunctionSpaceType &space, std::string name, typename DF::RangeType value )
    {
      return interpolant< DF >( space, std::move( name ), simpleGridFunction( space.gridPart(), [ value ] ( typename DF::DomainType ) { return value; }, 0 ) );
    }



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

      cls.def( "interpol", [] ( const Space &space, const GridFunction &gf, std::string name ) {
          return interpolant< DiscreteFunction >( space, name, gf );
        }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "interpol", [] ( const Space &space, typename Space::RangeType value, std::string name ) {
          return interpolant< DiscreteFunction >( space, name, value );
        }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "interpol", [] ( const Space &space, const GridFunction &gf ) {
          std::string name = gf.name();
          return interpolant< DiscreteFunction >( space, name, gf );
        }, pybind11::keep_alive< 0, 1 >() );
      cls.def( "interpol", [] ( const Space &space, typename Space::RangeType value ) {
          std::string name = "constant";
          return interpolant< DiscreteFunction >( space, name, value );
        }, pybind11::keep_alive< 0, 1 >() );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SPACE_HH
