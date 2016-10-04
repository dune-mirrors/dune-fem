#ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
#define DUNE_FEMPY_PY_GRID_GRIDPART_HH

#include <cassert>

#include <map>
#include <string>
#include <type_traits>
#include <utility>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/gridpart/common/gridpart2gridview.hh>

#include <dune/corepy/grid/hierarchical.hh>
#include <dune/corepy/grid/vtk.hh>

#include <dune/fempy/grid/gridpartadapter.hh>

#include <dune/corepy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {

    namespace detail
    {

      // GridModificationListener
      // ------------------------

      template< class Grid >
      class GridModificationListener final
        : public CorePy::GridModificationListener< Grid >
      {
        typedef Fem::DofManager< Grid > DofManager;

      public:
        GridModificationListener ( const Grid &grid )
          : dofManager_( DofManager::instance( grid ) )
        {
          std::cout << "created modification listener" << std::endl;
        }

        ~GridModificationListener ()
        {
          std::cout << "destroyed modification listener" << std::endl;
        }

        virtual void preModification ( const Grid &grid )
        {
          std::cout << "pre modification listener" << std::endl;
        }

        virtual void postModification ( const Grid &grid )
        {
          std::cout << "post modification listener" << std::endl;
          dofManager_.resize();
          dofManager_.compress();
        }

      private:
        DofManager &dofManager_;
      };


      template< class Grid >
      inline static void addGridModificationListener ( const Grid &grid )
      {
        typedef GridModificationListener< Grid > Listener;
        for( const auto &listener : CorePy::detail::gridModificationListeners( grid ) )
        {
          if( dynamic_cast< Listener * >( listener.second ) )
            return;
        }

        pybind11::handle nurse = pybind11::detail::get_object_handle( &grid, pybind11::detail::get_type_info( typeid( Grid ) ) );
        CorePy::detail::addGridModificationListener( grid, new Listener( grid ), nurse );
      }



      // GridPartConverter
      // -----------------

      template< class GV >
      struct GridPartConverter
      {
        typedef GV GridView;
        typedef GridPartAdapter< GV > GridPart;

        GridPart &operator() ( pybind11::handle gridView )
        {
          auto result = instances_.emplace( gridView.ptr(), nullptr );
          auto pos = result.first;
          if( result.second )
          {
            GridView view = gridView.template cast< GridView >();

            // add grid modification listener (if not registered)
            addGridModificationListener( view.grid() );

            // create new gridpart object
            pos->second = new GridPart( view );

            // create Python guard object, removing the grid part once the grid view dies
            pybind11::cpp_function remove_gridpart( [ this, pos ] ( pybind11::handle weakref ) {
                delete pos->second;
                instances_.erase( pos );
                weakref.dec_ref();
              } );
            pybind11::weakref weakref( gridView, remove_gridpart );
            weakref.release();
          }
          assert( pos->second );
          return *pos->second;
        }

      private:
        std::map< PyObject *, GridPart * > instances_;
      };


      template< class GP >
      struct GridPartConverter< Dune::GridView< Fem::GridPart2GridViewTraits< GP > > >
      {
        typedef GP GridPart;
        typedef Dune::GridView< Fem::GridPart2GridViewTraits< GP > > GridView;

        GridPart &operator() ( pybind11::handle gridView )
        {
          return const_cast< GridPart & >( gridView.template cast< GridView >().impl().gridPart() );
        }
      };


#if 0
      template< class GP >
      struct GridPartConverter< Fem::GridPart2GridView< GP > >
        : public GridPartConverter< Dune::GridView< Fem::GridPart2GridViewTraits< GP > > >
      {};
#endif



      // gridPartConverter
      // -----------------

      template< class GridView >
      inline GridPartConverter< GridView > &gridPartConverter ()
      {
        static GridPartConverter< GridView > converter;
        return converter;
      }

    } // namespace detail



    // GridPart
    // --------

    template< class GridView >
    using GridPart = typename detail::GridPartConverter< GridView >::GridPart;



    // gridPart
    // --------

    template< class GridView >
    inline static GridPart< GridView > &gridPart ( pybind11::handle gridView )
    {
      return detail::gridPartConverter< GridView >()( std::move( gridView ) );
    }



    // constructGridPart
    // -----------------

    template< class GridPart, class... Args >
    inline static void constructGridPart ( Dune::GridView< Fem::GridPart2GridViewTraits< GridPart > > &gridView,  Args &&... args )
    {
      typedef Dune::GridView< Fem::GridPart2GridViewTraits< GridPart > > GridView;

      GridPart *gridPart = new GridPart( std::forward< Args >( args )... );
      new (&gridView) GridView( static_cast< GridView >( *gridPart ) );

      // obtain Python object for grid view
      pybind11::handle nurse = pybind11::detail::get_object_handle( &gridView, pybind11::detail::get_type_info( typeid( GridView ) ) );
      if( !nurse )
        return;

      // create Python guard object, removing the grid part once the grid view dies
      pybind11::cpp_function remove_gridpart( [ gridPart ] ( pybind11::handle weakref ) {
          delete gridPart;
          weakref.dec_ref();
        } );
      pybind11::weakref weakref( nurse, remove_gridpart );
      weakref.release();
    }



#if 0
    // registerGridPart
    // ----------------

    namespace detail
    {

      template< class GridPart, class Field, int dimR >
      double l2Norm ( pybind11::object &o )
      {
        typedef Dune::FemPy::VirtualizedGridFunction< GridPart, Dune::FieldVector< Field, dimR > > VF;
        const VF &u = o.cast< VF >();
        Dune::Fem::L2Norm< GridPart > norm( u.gridPart(), 8 );
        return norm.norm( u );
      }

      template< class GridPart, int ... dimRange >
      auto defL2Norm( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange ... > )
      {
        typedef std::function< double (pybind11::object &) > Dispatch;
        std::array< Dispatch, sizeof ... ( dimRange ) > dispatch = {{ Dispatch( l2Norm< GridPart, double, dimRange > ) ... }};

        return [ dispatch ] ( pybind11::object gp, pybind11::object func ) {
           //const GridPart &gridPart = gp.cast< const GridPart & >();
           pybind11::object dimRobj = func.attr( "dimRange" );
           int dimR = func.attr( "dimRange" ).cast< int >();
           if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
             DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
           return dispatch[ static_cast< std::size_t >( dimR ) ]( func );
        };
      }

      template< class GridPart, class Cls >
      void registerGridPartConstructorFromGrid ( Cls &cls, std::false_type )
      {}

      template< class GridPart, class Cls >
      void registerGridPartConstructorFromGrid ( Cls &cls, std::true_type )
      {
        typedef typename GridPart::GridType Grid;
        cls.def( "__init__", [] ( GridPart &instance, Grid &grid ) {
            new (&instance)GridPart( grid );
          }, pybind11::keep_alive< 1, 2 >() );
      }

      template< class GridPart, class Cls >
      void registerGridPart ( pybind11::handle scope, Cls &cls )
      {
        typedef typename GridPart::GridType Grid;

        registerGridPartConstructorFromGrid< GridPart, Cls >( cls, std::is_constructible< GridPart, Grid & >() );
        cls.attr( "dimGrid" ) = pybind11::int_( GridPart::dimension );
        cls.attr( "dimWorld" ) = pybind11::int_( GridPart::dimensionworld );

        cls.def( "globalGridFunction", defGlobalGridFunction< GridPart >( cls, "GlobalGridFunction", std::make_integer_sequence< int, 11 >() ) );
        cls.def( "localGridFunction", defLocalGridFunction< GridPart >( cls, "LocalGridFunction", std::make_integer_sequence< int, 11 >() ) );

        cls.def( "l2Norm", defL2Norm< GridPart >( cls, "l2Norm", std::make_integer_sequence< int, 11 >() ) );
      }

    } // namespace detail
#endif // #if 0

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_GRIDPART_HH
