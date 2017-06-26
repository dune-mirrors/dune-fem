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

#include <dune/corepy/grid/numpy.hh>

#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // pointData
    // ---------

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, int level, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const int dimGrid = GridPart::dimension;

      const GeometryType gtSimplex( GeometryType::simplex, dimGrid );
      std::vector< Range > values;
      for( const auto &element : elements( static_cast< typename GridPart::GridViewType >( gridFunction.gridPart() ), ps ) )
      {
        const auto localFunction = gridFunction.localFunction( element );
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), gtSimplex );

        for( auto it = refinement.vBegin( level ), end = refinement.vEnd( level ); it != end; ++it )
        {
          Range value;
          localFunction.evaluate( it.coords(), value );
          values.push_back( value );
        }
      }

      return CorePy::makeNumPyArray< Field >( values, { values.size(), Range::dimension } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return pointData( gridFunction, 0, ps );
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    pointData ( const GridFunction &gridFunction, int level = 0 )
    {
      return pointData( gridFunction, level, Partitions::all );
    }



    // cellData
    // --------

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, int level, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const int dimGrid = GridPart::dimension;

      const GeometryType gtSimplex( GeometryType::simplex, dimGrid );
      std::vector< Range > values;
      for( const auto &element : elements( static_cast< typename GridPart::GridViewType >( gridFunction.gridPart() ), ps ) )
      {
        const auto localFunction = gridFunction.localFunction( element );
        const auto &refinement = buildRefinement< dimGrid, double >( element.type(), gtSimplex );

        for( auto it = refinement.eBegin( 0 ), end = refinement.eEnd( 0 ); it != end; ++it )
        {
          Range value;
          localFunction.evaluate( it.coords(), value );
          values.push_back( value );
        }
      }

      return CorePy::makeNumPyArray< Field >( values, { values.size(), Range::dimension } );
    }

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      return cellData( gridFunction, 0, ps );
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type >
    cellData ( const GridFunction &gridFunction, int level = 0 )
    {
      return cellData( gridFunction, level, Partitions::all );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
