#ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
#define DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>

#include <dune/common/stringutility.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/io/file/gmshreader.hh>

#include <dune/fempy/grid/adaptation.hh>
#include <dune/fempy/py/grid/entity.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>
#include <dune/fempy/py/common/common.hh>
#include <dune/fempy/pybind11/functional.h>
#include <dune/fempy/pybind11/numpy.h>
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/stl.h>
#include <dune/fempy/pybind11/extensions.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // gridAdaptationInstances
      // -----------------------

      template< class Grid >
      inline std::map< Grid *, std::unique_ptr< GridAdaptation< Grid > > > &gridAdaptationInstances ()
      {
        static std::map< Grid *, std::unique_ptr< GridAdaptation< Grid > > > instances;
        return instances;
      }



      // HierarchicalGridDeleter
      // -----------------------

      template< class Grid >
      struct HierarchicalGridDeleter
      {
        void operator() ( Grid *ptr )
        {
          auto &instances = gridAdaptationInstances< Grid >();
          auto it = instances.find( ptr );
          if( it != instances.end() )
            instances.erase( it );
          delete ptr;
        }
      };

    } // namespace detail



    // gridAdaptation
    // --------------

    template< class Grid >
    inline GridAdaptation< Grid > &gridAdaptation ( Grid &grid )
    {
      auto &instances = detail::gridAdaptationInstances< Grid >();
      std::unique_ptr< GridAdaptation< Grid > > &gridAdaptation = instances[ &grid ];
      if( !gridAdaptation )
        gridAdaptation.reset( new GridAdaptation< Grid >( grid ) );
      return *gridAdaptation;
    }

    template <class Grid>
    Grid *reader(const std::tuple<Reader,std::string> &constructor)
    {
      switch (std::get<0>(constructor))
      {
        case dgf:
        {
          GridPtr< Grid > gridPtr( std::get<1>(constructor) ) ;
          gridPtr->loadBalance();
          return gridPtr.release();
        }
        case dgfString:
        {
          std::istringstream in( std::get<1>(constructor) );
          GridPtr< Grid > gridPtr( in );
          gridPtr->loadBalance();
          return gridPtr.release();
        }
        case gmsh:
        {
#if 1
          std::cout << "PROBLEM HERE WITH NON STANDARD METHOD BEING USED IN GMSHREADER ON GRIDFACTORY - FAILS WITH YASPGRID"
            << std::endl;
#else
          Dune::GridFactory<Grid> gridFactory;
          Dune::GmshReader<Grid>::read(gridFactory, std::get<1>(constructor), false, false);
          return gridFactory.createGrid();
#endif
        }
      }
      return nullptr;
    }
    template <class Grid>
    Grid *reader(const pybind11::dict &constructor)
    {
      typedef typename Grid::ctype ctype;
      GridFactory< Grid > factory;

      // fill factory
      for (auto item : constructor)
      {
        std::string type = item.first.template cast<std::string>();
        pybind11::buffer_info buf = item.second.template cast< pybind11::buffer >().request();
        if (type == "vertex")
        {
          if( (buf.ndim != 2) || (buf.shape[ 1 ] != Grid::dimensionworld) )
            throw std::invalid_argument( "points array must be of shape (*, " + std::to_string( Grid::dimensionworld ) + ")" );
          for( std::size_t i = 0; i < buf.shape[ 0 ]; ++i )
          {
            const std::size_t offset = i * (buf.strides[ 0 ] / sizeof( ctype ));
            FieldVector< ctype, Grid::dimensionworld > x;
            for( int j = 0; j < Grid::dimensionworld; ++j )
              x[ j ] = static_cast< ctype >( static_cast< ctype * >( buf.ptr )[ offset + j * (buf.strides[ 1 ] / sizeof( ctype )) ] );
            factory.insertVertex( x );
          }
        }
        else if (type == "simplex")
        {
          if( (buf.ndim != 2) || (buf.shape[ 1 ] != Grid::dimension+1) )
            throw std::invalid_argument( "simplices array must be of shape (*, " + std::to_string( Grid::dimension+1 ) + ")" );
          GeometryType type( GeometryType::simplex, Grid::dimension );
          std::vector< unsigned int > vertices( Grid::dimension+1 );
          for( std::size_t i = 0; i < buf.shape[ 0 ]; ++i )
          {
            const std::size_t offset = i * (buf.strides[ 0 ] / sizeof( int ));
            for( int j = 0; j <= Grid::dimension; ++j )
              vertices[ j ] = static_cast< int * >( buf.ptr )[ offset + j * (buf.strides[ 1 ] / sizeof( int )) ];
            factory.insertElement( type, vertices );
          }
        }
        else if (type == "cube")
        {
          unsigned int size = 2 << Grid::dimension;
          if( (buf.ndim != 2) || (buf.shape[ 1 ] != size) )
            throw std::invalid_argument( "simplices array must be of shape (*, " + std::to_string( size ) + ")" );
          GeometryType type( GeometryType::cube, Grid::dimension );
          std::vector< unsigned int > vertices( pow(2,Grid::dimension) );
          for( std::size_t i = 0; i < buf.shape[ 0 ]; ++i )
          {
            const std::size_t offset = i * (buf.strides[ 0 ] / sizeof( int ));
            for( int j = 0; j <= Grid::dimension; ++j )
              vertices[ j ] = static_cast< int * >( buf.ptr )[ offset + j * (buf.strides[ 1 ] / sizeof( int )) ];
            factory.insertElement( type, vertices );
          }
        }
      }
      // create grid
      return factory.createGrid();
    }

    // registerHierarchicalGrid
    // ------------------------

    template< class Grid >
    inline auto registerHierarchicalGrid ( pybind11::handle scope )
    {
      typedef GridAdaptation< Grid > Adaptation;

      registerGridEntities< Grid >( scope );

      typedef typename Adaptation::Element Element;
      typedef typename Adaptation::Marker Marker;

      typedef VirtualizedRestrictProlong< Grid > RestrictProlong;

      typedef detail::HierarchicalGridDeleter< Grid > Deleter;
      pybind11::class_< Grid, std::unique_ptr< Grid, Deleter > > cls( scope, "HierarchicalGrid" );

      cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( Grid::dimension ) );
      cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( Grid::dimensionworld ) );

      cls.def( "__repr__", [] ( const Grid &grid ) -> std::string { return "HierarchicalGrid"; } );

      cls.def( "globalRefine", [] ( Grid &grid ) { gridAdaptation( grid ).globalRefine( 1 ); } );
      cls.def( "globalRefine", [] ( Grid &grid, int level ) { gridAdaptation( grid ).globalRefine( level ); } );

      detail::clsVirtualizedRestrictProlong< Grid >( cls );

      pybind11::enum_< Marker > marker( cls, "marker" );
      marker.value( "coarsen", Marker::Coarsen );
      marker.value( "keep", Marker::Keep );
      marker.value( "refine", Marker::Refine );

      cls.def( "mark", [] ( Grid &grid, const std::function< Marker( const Element &e ) > &marker ) {
          return gridAdaptation( grid ).mark( marker );
        } );

      cls.def( "adapt", [] ( Grid &grid ) {
          std::array< RestrictProlong, 0 > rpList;
          gridAdaptation( grid ).adapt( rpList.begin(), rpList.end() );
        } );


      cls.def( "adapt", [] ( Grid &grid, const std::list< RestrictProlong > &rpList ) {
          std::cout << "adapting grid and " << rpList.size() << " functions..." << std::endl;
          gridAdaptation( grid ).adapt( rpList.begin(), rpList.end() );
        } );

      cls.def( "loadBalance", [] ( Grid &grid ) {
          std::array< RestrictProlong, 0 > rpList;
          gridAdaptation( grid ).loadBalance( rpList.begin(), rpList.end() );
        } );

      cls.def( "loadBalance", [] ( Grid &grid, const std::list< RestrictProlong > &rpList ) {
          std::cout << "loadbalanding grid and " << rpList.size() << " functions..." << std::endl;
          gridAdaptation( grid ).loadBalance( rpList.begin(), rpList.end() );
        } );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
