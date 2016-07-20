#ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
#define DUNE_FEMPY_PY_GRID_NUMPY_HH

#include <dune/common/ftraits.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/quadrature/cornerpointset.hh>

#include <dune/fempy/pybind11/numpy.h>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // coordinates
    // -----------

    template< class GridView, class Mapper >
    inline static pybind11::array_t< typename GridView::ctype > coordinates ( const GridView &gridView, const Mapper &mapper )
    {
      typedef typename GridView::ctype ctype;

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( mapper.size() ), static_cast< std::size_t >( GridView::dimensionworld ) };
      const std::vector< std::size_t > stride{ GridView::dimensionworld * sizeof( ctype ), sizeof( ctype ) };
      pybind11::array_t< ctype > coords( pybind11::buffer_info( nullptr, sizeof( ctype ), pybind11::format_descriptor< ctype >::value, 2, shape, stride ) );

      pybind11::buffer_info info = coords.request( true );
      for( const auto &vertex : vertices( gridView, Partitions::all ) )
      {
        typename Mapper::Index index;
        if( !mapper.contains( vertex, index ) )
          continue;

        const auto x = vertex.geometry().center();
        std::copy( x.begin(), x.end(), static_cast< ctype * >( info.ptr ) + GridView::dimensionworld * index );
      }

      return coords;
    }

    template< class GridView >
    inline static pybind11::array_t< typename GridView::ctype > coordinates ( const GridView &gridView )
    {
      MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout > mapper( gridView );
      return coordinates( gridView, mapper );
    }



    // tessellate
    // ----------

    template< class GridView, class Mapper, unsigned int partitions >
    inline static pybind11::array_t< int > tesselate ( const GridView &gridView, const Mapper &mapper, PartitionSet< partitions > ps )
    {
      const int dimension = GridView::dimension;

      std::vector< std::array< int, dimension+1 > > simplices;
      for( const auto &element : elements( gridView, ps ) )
      {
        if( !element.type().isSimplex() )
          DUNE_THROW( NotImplemented, "Only simplicial grids can be tessealted right now." );

        std::array< int, dimension+1 > simplex;
        for( int i = 0; i <= dimension; ++i )
          simplex[ i ] = mapper.subIndex( element, i, dimension );

        simplices.push_back( simplex );
      }

      const std::vector< std::size_t > shape{ simplices.size(), static_cast< std::size_t >( dimension+1 ) };
      const std::vector< std::size_t > stride{ (dimension+1) * sizeof( int ), sizeof( int ) };
      pybind11::array_t< int > result( pybind11::buffer_info( nullptr, sizeof( int ), pybind11::format_descriptor< int >::value, 2, shape, stride ) );

      pybind11::buffer_info info = result.request( true );
      int *out = static_cast< int * >( info.ptr );
      for( const auto &simplex : simplices )
        out = std::copy( simplex.begin(), simplex.end(), out );

      return result;
    }

    template< class GridView, unsigned int partitions >
    inline static pybind11::array_t< int > tesselate ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout > mapper( gridView );
      return tesselate( gridView, ps );
    }

    template< class GridView >
    inline static pybind11::array_t< int > tesselate ( const GridView &gridView )
    {
      return tesselate( gridView, Partitions::all );
    }



    // pointData
    // ---------

    template< class GridFunction, class Mapper >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > pointData ( const GridFunction &gridFunction, const Mapper &mapper )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( mapper.size() ), static_cast< std::size_t >( Range::dimension ) };
      const std::vector< std::size_t > stride{ Range::dimension * sizeof( Field ), sizeof( Field ) };
      pybind11::array_t< Field > data( pybind11::buffer_info( nullptr, sizeof( Field ), pybind11::format_descriptor< Field >::value, 2, shape, stride ) );

      pybind11::buffer_info info = data.request( true );
      for( const auto &element : elements( static_cast< typename GridPart::GridViewType >( gridFunction.gridPart() ), Partitions::all ) )
      {
        const auto localFunction = gridFunction.localFunction( element );

        Fem::CornerPointSet< GridPart > corners( element );
        for( std::size_t i = 0; i < corners.nop(); ++i )
        {
          typename Mapper::Index index;
          if( !mapper.contains( element, i, GridPart::dimension, index ) )
            continue;
          Range value;
          localFunction.evaluate( corners[ i ], value );
          std::copy( value.begin(), value.end(), static_cast< Field * >( info.ptr ) + Range::dimension * index );
        }
      }

      return data;
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > pointData ( const GridFunction &gridFunction )
    {
      typedef typename GridFunction::GridPartType::GridViewType GridView;
      MultipleCodimMultipleGeomTypeMapper< GridView, MCMGVertexLayout > mapper( static_cast< GridView >( gridFunction.gridPart() ) );
      return pointData( gridFunction, mapper );
    }



    // cellData
    // --------

    template< class GridFunction, unsigned int partitions >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > cellData ( const GridFunction &gridFunction, PartitionSet< partitions > ps )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::LocalFunctionType LocalFunction;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      std::vector< Range > cellData;
      Fem::LocalAverage< LocalFunction, GridPart > localAverage;
      for( const auto &element : elements( static_cast< typename GridPart::GridViewType >( gridFunction.gridPart() ), ps ) )
      {
        if( !element.type().isSimplex() )
          DUNE_THROW( NotImplemented, "Only simplicial grids can be tessealted right now." );

        Range value;
        localAverage( gridFunction.localFunction( element ), value );
        cellData.push_back( value );
      }

      const std::vector< std::size_t > shape{ cellData.size(), static_cast< std::size_t >( Range::dimension ) };
      const std::vector< std::size_t > stride{ Range::dimension * sizeof( Field ), sizeof( Field ) };
      pybind11::array_t< Field > result( pybind11::buffer_info( nullptr, sizeof( Field ), pybind11::format_descriptor< Field >::value, 2, shape, stride ) );

      pybind11::buffer_info info = result.request( true );
      Field *out = static_cast< Field * >( info.ptr );
      for( const Range &value : cellData )
        out = std::copy( value.begin(), value.end(), out );

      return result;
    }

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > cellData ( const GridFunction &gridFunction )
    {
      return cellData( gridFunction, Partitions::all );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
