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

#include <dune/corepy/grid/hierarchical.hh>

#include <dune/fempy/grid/adaptation.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>

#include <dune/corepy/pybind11/functional.h>
#include <dune/corepy/pybind11/numpy.h>
#include <dune/corepy/pybind11/pybind11.h>
#include <dune/corepy/pybind11/stl.h>
#include <dune/corepy/pybind11/extensions.h>

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

    // registerHierarchicalGrid
    // ------------------------

    template< class Grid >
    inline auto registerHierarchicalGrid ( pybind11::module scope )
    {
      typedef detail::HierarchicalGridDeleter< Grid > Deleter;
      auto cls = Dune::CorePy::registerHierarchicalGrid<Grid, std::unique_ptr<Grid,Deleter> >(scope);

#if 0
      pybind11::class_< Grid, std::unique_ptr< Grid, Deleter > > cls( scope, "HierarchicalGrid" );

      cls.attr( "dimGrid" ) = pybind11::int_( static_cast< int >( Grid::dimension ) );
      cls.attr( "dimWorld" ) = pybind11::int_( static_cast< int >( Grid::dimensionworld ) );
      cls.def( "__repr__", [] ( const Grid &grid ) -> std::string { return "HierarchicalGrid"; } );

      module.def( "reader", [](const std::tuple<Reader,std::string> &constructor)
          { return reader<Grid>(constructor); } );
      module.def( "reader", [](const std::string &constructor)
          { return reader<Grid>(std::make_tuple(Reader::dgf,constructor)); } );
      module.def( "reader", [](const pybind11::dict &constructor)
          { return reader<Grid>(constructor); } );

      registerGridEntities< Grid >( scope );
#endif

      // we provide extra global refine methods for fem - in grid.py these
      // are called in the globalRefine methods exported on the leafGrid
      cls.def( "femGlobalRefine", [] ( Grid &grid ) { gridAdaptation( grid ).globalRefine( 1 ); } );
      cls.def( "femGlobalRefine", [] ( Grid &grid, int level ) { std::cout << "RIGHT GLOBALREFINE" << std::endl; gridAdaptation( grid ).globalRefine( level ); } );

      detail::clsVirtualizedRestrictProlong< Grid >( cls );

      typedef VirtualizedRestrictProlong< Grid > RestrictProlong;
      typedef GridAdaptation< Grid > Adaptation;
      typedef typename Adaptation::Marker Marker;
      typedef typename Adaptation::Element Element;

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
