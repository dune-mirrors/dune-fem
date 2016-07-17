#ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
#define DUNE_FEMPY_PY_GRID_NUMPY_HH

#include <dune/common/ftraits.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dune/fem/function/localfunction/average.hh>
#include <dune/fem/quadrature/cornerpointset.hh>

#include <dune/corepy/pybind11/numpy.h>
#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    // coordinates
    // -----------

    template< class GridView >
    inline static pybind11::array_t< typename GridView::ctype > coordinates ( const GridView &gridView )
    {
      typedef typename GridView::ctype ctype;
      const auto &indexSet = gridView.indexSet();

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( indexSet.size( GridView::dimension ) ), static_cast< std::size_t >( GridView::dimensionworld ) };
      const std::vector< std::size_t > stride{ GridView::dimensionworld * sizeof( ctype ), sizeof( ctype ) };
      pybind11::array_t< ctype > coords( pybind11::buffer_info( nullptr, sizeof( ctype ), pybind11::format_descriptor< ctype >::value, 2, shape, stride ) );

      pybind11::buffer_info info = coords.request( true );
      for( const auto &vertex : vertices( gridView, Partitions::all ) )
      {
        const int index = indexSet.index( vertex );
        const auto x = vertex.geometry().center();
        std::copy( x.begin(), x.end(), static_cast< ctype * >( info.ptr ) + GridView::dimensionworld * index );
      }

      return coords;
    }



    // tessellate
    // ----------

    // Warning: Only Partitions::all is working right now.
    template< class GridView, unsigned int partitions >
    inline static pybind11::array_t< int > tesselate ( const GridView &gridView, PartitionSet< partitions > ps )
    {
      const int dimension = GridView::dimension;

      const auto &indexSet = gridView.indexSet();

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( indexSet.size( 0 ) ), static_cast< std::size_t >( dimension+1 ) };
      const std::vector< std::size_t > stride{ (dimension+1) * sizeof( int ), sizeof( int ) };
      pybind11::array_t< int > simplices( pybind11::buffer_info( nullptr, sizeof( int ), pybind11::format_descriptor< int >::value, 2, shape, stride ) );

      pybind11::buffer_info info = simplices.request( true );
      for( const auto &element : elements( gridView, ps ) )
      {
        const int index = indexSet.index( element );
        if( !element.type().isSimplex() )
          DUNE_THROW( NotImplemented, "Only simplicial grids can be tessealted right now." );
        for( int i = 0; i <= dimension; ++i )
          static_cast< int * >( info.ptr )[ (dimension+1)*index + i ] = indexSet.subIndex( element, i, dimension );
      }

      return simplices;
    }

    template< class GridView >
    inline static pybind11::array_t< int > tesselate ( const GridView &gridView )
    {
      return tesselate( gridView, Partitions::all );
    }



    // pointData
    // ---------

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > pointData ( const GridFunction &gridFunction )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      const auto &indexSet = gridFunction.gridPart().indexSet();

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( indexSet.size( GridPart::dimension ) ), static_cast< std::size_t >( Range::dimension ) };
      const std::vector< std::size_t > stride{ Range::dimension * sizeof( Field ), sizeof( Field ) };
      pybind11::array_t< Field > data( pybind11::buffer_info( nullptr, sizeof( Field ), pybind11::format_descriptor< Field >::value, 2, shape, stride ) );

      pybind11::buffer_info info = data.request( true );
      for( const auto &element : elements( static_cast< typename GridPart::GridViewType >( gridFunction.gridPart() ), Partitions::all ) )
      {
        const auto localFunction = gridFunction.localFunction( element );

        Fem::CornerPointSet< GridPart > corners( element );
        for( std::size_t i = 0; i < corners.nop(); ++i )
        {
          const std::size_t index = indexSet.subIndex( element, i, GridPart::dimension );
          Range value;
          localFunction.evaluate( corners[ i ], value );
          std::copy( value.begin(), value.end(), static_cast< Field * >( info.ptr ) + Range::dimension * index );
        }
      }

      return data;
    }



    // cellData
    // --------

    template< class GridFunction >
    inline static pybind11::array_t< typename FieldTraits< typename GridFunction::RangeType >::field_type > cellData ( const GridFunction &gridFunction )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::LocalFunctionType LocalFunction;
      typedef typename GridFunction::RangeType Range;
      typedef typename FieldTraits< Range >::field_type Field;

      typedef typename GridPart::GridViewType GridView;
      MultipleCodimMultipleGeomTypeMapper< GridView, MCMGElementLayout > mapper( static_cast< GridView >( gridFunction.gridPart() ) );

      const std::vector< std::size_t > shape{ static_cast< std::size_t >( mapper.size() ), static_cast< std::size_t >( Range::dimension ) };
      const std::vector< std::size_t > stride{ Range::dimension * sizeof( Field ), sizeof( Field ) };
      pybind11::array_t< Field > data( pybind11::buffer_info( nullptr, sizeof( Field ), pybind11::format_descriptor< Field >::value, 2, shape, stride ) );

      pybind11::buffer_info info = data.request( true );
      Fem::LocalAverage< LocalFunction, GridPart > localAverage;
      for( const auto &element : elements( static_cast< GridView >( gridFunction.gridPart() ), Partitions::all ) )
      {
        Range value;
        localAverage( gridFunction.localFunction( element ), value );
        const int index = mapper.index( element );
        std::copy( value.begin(), value.end(), static_cast< Field * >( info.ptr ) + Range::dimension * index );
      }

      return data;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_NUMPY_HH
