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
      auto cls = pybind11::class_<Grid, std::unique_ptr<Grid,Deleter>>(scope, "HierarchicGrid");
      // we provide extra global refine methods for fem - in grid.py these
      // are called in the globalRefine methods exported on the leafGrid
      cls.def( "globalRefine", [] ( Grid &grid ) { gridAdaptation( grid ).globalRefine( 1 ); } );
      cls.def( "globalRefine", [] ( Grid &grid, int level ) { gridAdaptation( grid ).globalRefine( level ); } );

      detail::clsVirtualizedRestrictProlong< Grid >( cls );

      typedef VirtualizedRestrictProlong< Grid > RestrictProlong;
      typedef GridAdaptation< Grid > Adaptation;

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

      Dune::CorePy::registerHierarchicalGrid<Grid>(scope,cls);

      return cls;
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_HIERARCHICAL_HH
