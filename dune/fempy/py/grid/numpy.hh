#ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
#define DUNE_FEMPY_PY_GRID_NUMPY_HH

#include <cstddef>

#include <algorithm>
#include <vector>

#include <dune/common/ftraits.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/virtualrefinement.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/python/grid/numpy.hh>

#include <dune/fem/function/localfunction/const.hh>

#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // pointData
    // ---------

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, RefinementIntervals intervals, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const int dimGrid = GridPart::dimension;

      std::vector< Range > values;
      Dune::Fem::ConstLocalFunction<GridFunction> localFunction(gridFunction);
      for( const auto &element : elements( gridFunction.gridPart(), ps ) )
      {
        localFunction.bind(element);
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), GeometryTypes::simplex( dimGrid ) );

        for( auto it = refinement.vBegin( intervals ), end = refinement.vEnd( intervals ); it != end; ++it )
        {
          Range value;
          localFunction.evaluate( it.coords(), value );
          values.push_back( value );
        }
      }

      return Python::makeNumPyArray< Field >( values, { values.size(), Range::dimension } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return pointData( gridFunction, refinementLevels( 0 ), ps );
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, RefinementIntervals intervals = refinementLevels( 0 ) )
    {
      return pointData( gridFunction, intervals, Partitions::all );
    }



    // cellData
    // --------

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, RefinementIntervals intervals, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const int dimGrid = GridPart::dimension;

      std::vector< Range > values;
      Dune::Fem::ConstLocalFunction<GridFunction> localFunction(gridFunction);
      for( const auto &element : elements( gridFunction.gridPart(), ps ) )
      {
        localFunction.bind( element );
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), GeometryTypes::simplex( dimGrid ) );

        for( auto it = refinement.eBegin( intervals ), end = refinement.eEnd( intervals ); it != end; ++it )
        {
          Range value;
          localFunction.evaluate( it.coords(), value );
          values.push_back( value );
        }
      }

      return Python::makeNumPyArray< Field >( values, { values.size(), Range::dimension } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return cellData( gridFunction, refinementIntervals( 0 ), ps );
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, RefinementIntervals intervals = refinementLevels( 0 ) )
    {
      return cellData( gridFunction, intervals, Partitions::all );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
