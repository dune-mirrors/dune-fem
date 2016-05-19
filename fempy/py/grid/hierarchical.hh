#ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
#define DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH

#include <array>
#include <functional>
#include <list>
#include <map>
#include <memory>

#include <dune/fempy/grid/hierarchical.hh>
#include <dune/fempy/py/grid/entity.hh>
#include <dune/fempy/pybind11/functional.h>
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/stl.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // hierarchicalGridInstances
      // -------------------------

      template< class Grid >
      inline std::map< Grid *, HierarchicalGrid< Grid > > &hierarchicalGridInstances ()
      {
        static std::map< Grid *, HierarchicalGrid< Grid > > instances;
        return instances;
      }



      // HierarchicalGridDeleter
      // -----------------------

      template< class HierarchicalGrid >
      struct HierarchicalGridDeleter
      {
        void operator() ( HierarchicalGrid *ptr )
        {
          auto &instances = hierarchicalGridInstances< typename HierarchicalGrid::Grid >();
          auto it = instances.find( ptr->grid().get() );
          if( it != instances.end() )
            instances.erase( it );
          delete ptr;
        }
      };

    } // namespace detail



    // hierarchicalGrid
    // ----------------

    template< class Grid >
    inline HierarchicalGrid< Grid > hierarchicalGrid ( Grid &grid )
    {
      const auto &instances = detail::hierarchicalGridInstances< Grid >();
      auto it = instances.find( &grid );
      if( it == instances.end() )
        throw std::invalid_argument( "Unknown hierarchical grid" );
      return it->second;
    }



    // registerHierarchicalGrid
    // ------------------------

    template< class HierarchicalGrid >
    inline auto registerHierarchicalGrid ( pybind11::handle scope )
    {
      typedef typename HierarchicalGrid::Grid Grid;

      registerGridEntities< Grid >( scope );

      typedef typename HierarchicalGrid::Element Element;
      typedef typename HierarchicalGrid::Marker Marker;

      typedef AdaptiveDofVector< Grid, double > GridFunction;

      typedef detail::HierarchicalGridDeleter< HierarchicalGrid > Deleter;
      pybind11::class_< HierarchicalGrid, std::unique_ptr< HierarchicalGrid, Deleter > > cls( scope, "HierarchicalGrid" );
      cls.def( "__init__", [] ( HierarchicalGrid &instance, std::string dgf ) {
          new (&instance) HierarchicalGrid( dgf );
          detail::hierarchicalGridInstances< Grid >().insert( std::make_pair( instance.grid().get(), instance ) );
        } );
      cls.def( "__repr__", [] ( const HierarchicalGrid &grid ) -> std::string { return "HierarchicalGrid"; } );

      cls.def( "globalRefine", &HierarchicalGrid::globalRefine );

      pybind11::enum_< Marker > marker( cls, "marker" );
      marker.value( "coarsen", HierarchicalGrid::Marker::Coarsen );
      marker.value( "keep", HierarchicalGrid::Marker::Keep );
      marker.value( "refine", HierarchicalGrid::Marker::Refine );

      cls.def( "mark", [] ( HierarchicalGrid &grid, const std::function< Marker( const Element &e ) > &marker ) {
          grid.mark( marker );
        } );

      cls.def( "adapt", [] ( HierarchicalGrid &grid ) {
          std::array< std::shared_ptr< GridFunction >, 0 > dfList;
          grid.adapt( dfList.begin(), dfList.end() );
        } );


      cls.def( "adapt", [] ( HierarchicalGrid &grid, const std::list< std::shared_ptr< GridFunction > > &dfList ) {
          std::cout << "adapting grid and " << dfList.size() << " functions..." << std::endl;
          grid.adapt( dfList.begin(), dfList.end() );
        } );

      cls.def( "loadBalance", [] ( HierarchicalGrid &grid ) {
          std::array< std::shared_ptr< GridFunction >, 0 > dfList;
          grid.loadBalance( dfList.begin(), dfList.end() );
        } );

      cls.def( "loadBalance", [] ( HierarchicalGrid &grid, const std::list< std::shared_ptr< GridFunction > > &dfList ) {
          std::cout << "loadbalanding grid and " << dfList.size() << " functions..." << std::endl;
          grid.loadBalance( dfList.begin(), dfList.end() );
        } );

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
