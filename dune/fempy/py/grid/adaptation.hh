#ifndef DUNE_FEMPY_PY_GRID_ADAPTATION_HH
#define DUNE_FEMPY_PY_GRID_ADAPTATION_HH

#include <cassert>

#include <array>
#include <list>
#include <map>

#include <dune/fem/storage/singleton.hh>

#include <dune/fempy/grid/adaptation.hh>
#include <dune/fempy/py/grid/restrictprolong.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // gridAdaptationInstances
      // -----------------------

      template< class Grid >
      inline std::map< Grid *, GridAdaptation< Grid > * > &gridAdaptationInstances ()
      {
        typedef std::map< Grid *, GridAdaptation< Grid > * > InstancesMapType;
        return Dune::Fem::Singleton< InstancesMapType >::instance();
      }

    } // namespace detail



    // gridAdaptation
    // --------------

    template< class Grid >
    inline static GridAdaptation< Grid > &gridAdaptation ( Grid &grid )
    {
      auto result = detail::gridAdaptationInstances< Grid >().insert( std::make_pair( &grid, nullptr ) );
      auto pos = result.first;
      if( result.second )
      {
        pos->second = new GridAdaptation< Grid >( grid );

        pybind11::handle nurse = pybind11::detail::get_object_handle( &grid, pybind11::detail::get_type_info( typeid( Grid ) ) );
        if( nurse )
        {
          pybind11::cpp_function remove_adaptation( [ pos ] ( pybind11::handle weakref ) {
              delete pos->second;
              detail::gridAdaptationInstances< Grid >().erase( pos );
              weakref.dec_ref();
            } );
          pybind11::weakref( nurse, remove_adaptation ).release();
        }
      }
      assert( pos->second );
      return *pos->second;
    }



    // registerGridAdaptation
    // ----------------------

    template< class Grid, class... options >
    inline static void registerGridAdaptation ( pybind11::module module, pybind11::class_< GridAdaptation< Grid >, options... > cls )
    {
      detail::clsVirtualizedRestrictProlong< Grid >( cls );

      typedef VirtualizedRestrictProlong< Grid > RestrictProlong;

      cls.def( "globalRefine", [] ( GridAdaptation< Grid > &self, int level, const std::list< RestrictProlong > &rpList ) {
          for (int i=0;i<level;++i)
          {
            self.markAll();
            self.adapt( rpList.begin(), rpList.end() );
          }
        } );

      cls.def( "adapt", [] ( GridAdaptation< Grid > &self ) {
          std::array< RestrictProlong, 0 > rpList;
          self.adapt( rpList.begin(), rpList.end() );
        } );

      cls.def( "adapt", [] ( GridAdaptation< Grid > &self, const std::list< RestrictProlong > &rpList ) {
          self.adapt( rpList.begin(), rpList.end() );
        } );

      cls.def( "loadBalance", [] ( GridAdaptation< Grid > &self ) {
          std::array< RestrictProlong, 0 > rpList;
          self.loadBalance( rpList.begin(), rpList.end() );
        } );

      cls.def( "loadBalance", [] ( GridAdaptation< Grid > &self, const std::list< RestrictProlong > &rpList ) {
          // if( self.grid().comm().rank() == 0 )
          //   std::cout << "loadbalanding grid and " << rpList.size() << " functions..." << std::endl;
          self.loadBalance( rpList.begin(), rpList.end() );
        } );

      module.def( "gridAdaptation", &gridAdaptation< Grid >, pybind11::return_value_policy::reference );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_ADAPTATION_HH
