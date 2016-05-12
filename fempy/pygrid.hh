#ifndef DUNE_FEMPY_PYGRID_HH
#define DUNE_FEMPY_PYGRID_HH

#include <list>
#include <memory>
#include <tuple>

#include <dune/common/std/utility.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/grid.hh>
#include <dune/fempy/pybind11/functional.h>
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/stl.h>
#include <dune/fempy/pygridfunction.hh>
#include <dune/fempy/pyvtk.hh>
#include <dune/fempy/vtk.hh>

namespace Dune
{

  namespace FemPy
  {

#if 0
    // registerIteratorRange
    // ---------------------

    template< class IteratorRange >
    void registerIteratorRange ( pybind11::handle scope, const char *name )
    {
      pybind11::class_< IteratorRange > cls( scope, name );
      cls.def( "__iter__", [] ( const IteratorRange &rg ) { return pybind11::make_iterator( rg.begin(), rg.end() ); }, pybind11::keep_alive< 0, 1 >() );
    }
#endif



    // PyCornerRange
    // -------------

    template< class Geometry >
    struct PyCorners
    {
      PyCorners ( const Geometry &geometry, pybind11::object ref )
        : geometry_( geometry ), ref_( ref )
      {}

      const Geometry &geometry () { return geometry_; }

    private:
      const Geometry &geometry_;
      pybind11::object ref_;
    };


    // PyCornerIterator
    // ----------------

    template< class Geometry >
    struct PyCornerIterator
    {
      PyCornerIterator ( const PyCorners< Geometry > corners ) : corners_( corners ) {}

      typename Geometry::GlobalCoordinate next ()
      {
        if( index_ == corners_.geometry().corners() )
          throw pybind11::stop_iteration();
        return corners_.geometry().corner( index_++ );
      }

    private:
      PyCorners< Geometry > corners_;
      int index_ = 0;
    };



    // registerGridGeometry
    // --------------------

    template< class Geometry >
    void registerGridGeometry ( pybind11::module module )
    {
      static const std::string clsName = "Geometry" + std::to_string( Geometry::mydimension );
      pybind11::class_< Geometry > cls( module, clsName.c_str() );

      pybind11::class_< PyCornerIterator< Geometry > > itCls( cls, "CornerIterator" );
      itCls.def( "__iter__", [] ( PyCornerIterator< Geometry > &it ) -> PyCornerIterator< Geometry > & { return it; } );
      itCls.def( "__next__", &PyCornerIterator< Geometry >::next );

      pybind11::class_< PyCorners< Geometry > > cCls( cls, "Corners" );
      cCls.def( "__iter__", [] ( const PyCorners< Geometry > &c ) { return PyCornerIterator< Geometry >( c ); } );

      cls.def_property_readonly( "corners", [] ( pybind11::object geo ) {
          return PyCorners< Geometry >( geo.cast< const Geometry & >(), geo );
        } );

      cls.def_property_readonly( "center", &Geometry::center );
      cls.def_property_readonly( "volume", &Geometry::volume );

      cls.def_property_readonly( "affine", &Geometry::affine );

      cls.def( "position", &Geometry::global );
      cls.def( "integrationElement", &Geometry::integrationElement );
    }


    template< class Entity >
    void registerGridEntity ( pybind11::module module )
    {
      registerGridGeometry< typename Entity::Geometry >( module );

      static const std::string entityName = "Entity" + std::to_string( Entity::codimension );
      pybind11::class_< Entity > entity( module, entityName.c_str() );

      entity.def_property_readonly_static( "codimension", [] () { return  Entity::codimension; } );

      entity.def_property_readonly( "geometry", &Entity::geometry );
      entity.def_property_readonly( "level", &Entity::level );
    }


    template< class Grid, int... codim >
    void registerGridEntities ( pybind11::module module, Std::integer_sequence< int, codim... > )
    {
      std::ignore = std::make_tuple( (registerGridEntity< typename Grid::template Codim< codim >::Entity >( module ), codim)... );
    }



    // registerHierarchicalGrid
    // ------------------------

    template< class Grid >
    void registerHierarchicalGrid ( pybind11::module module )
    {
      static auto common = pybind11::module::import( "dune.common" );
      static auto femmpi = pybind11::module::import( "dune.femmpi" );

      registerGridEntities< Grid >( module, Std::make_integer_sequence< int, Grid::dimension+1 >() );

      typedef HierarchicalGrid< Grid > HG;
      typedef typename HG::Element Element;
      typedef typename HG::Marker Marker;

      typedef AdaptiveGridFunction< Grid > GridFunction;

      pybind11::class_< HG > hg( module, "HierarchicalGrid" );
      hg.def( pybind11::init< std::string >() );
      hg.def( "__repr__", [] ( const HG &grid ) -> std::string { return "HierarchicalGrid"; } );

      hg.def( "globalRefine", &HG::globalRefine );

      pybind11::enum_< Marker > marker( hg, "marker" );
      marker.value( "coarsen", HG::Marker::Coarsen );
      marker.value( "keep", HG::Marker::Keep );
      marker.value( "refine", HG::Marker::Refine );

      hg.def( "mark", [] ( HG &grid, const std::function< Marker( const Element &e ) > &marker ) {
          grid.mark( marker );
        } );

      hg.def( "adapt", [] ( HG &grid ) {
          std::array< std::shared_ptr< GridFunction >, 0 > dfList;
          grid.adapt( dfList.begin(), dfList.end() );
        } );


      hg.def( "adapt", [] ( HG &grid, const std::list< std::shared_ptr< GridFunction > > &dfList ) {
          std::cout << "adapting grid and " << dfList.size() << " functions..." << std::endl;
          grid.adapt( dfList.begin(), dfList.end() );
        } );

      hg.def( "loadBalance", [] ( HG &grid ) {
          std::array< std::shared_ptr< GridFunction >, 0 > dfList;
          grid.loadBalance( dfList.begin(), dfList.end() );
        } );

      hg.def( "loadBalance", [] ( HG &grid, const std::list< std::shared_ptr< GridFunction > > &dfList ) {
          std::cout << "loadbalanding grid and " << dfList.size() << " functions..." << std::endl;
          grid.loadBalance( dfList.begin(), dfList.end() );
        } );
    }



    // PyGridPartRange
    // ---------------

    template< class GridPart, int codim >
    struct PyGridPartRange
    {
      typedef typename GridPart::template Codim< codim >::IteratorType Iterator;

      PyGridPartRange ( const GridPart &gridPart, pybind11::object ref )
        : gridPart_( gridPart ), ref_( std::move( ref ) )
      {}

      Iterator begin () const { return gridPart_.template begin< codim >(); }
      Iterator end () const { return gridPart_.template end< codim >(); }

    private:
      const GridPart &gridPart_;
      pybind11::object ref_;
    };



    // PyGridPartIterator
    // ------------------

    template< class GridPart, int codim >
    struct PyGridPartIterator
    {
      typedef PyGridPartRange< GridPart, codim > Range;
      typedef typename GridPart::template Codim< codim >::EntityType Entity;

      PyGridPartIterator ( const Range &range ) : range_( range ), it_( range_.begin() ) {}

      Entity next ()
      {
        if( it_ == range_.end() )
          throw pybind11::stop_iteration();

        Entity entity = *it_;
        ++it_;
        return entity;
      }

    private:
      Range range_;
      typename Range::Iterator it_;
    };



    // registerPyGridPartRange
    // -----------------------

    template< class GridPart, int codim >
    void registerPyGridPartRange ( pybind11::handle scope, const char *rgName )
    {
      typedef PyGridPartRange< GridPart, codim > Range;
      typedef PyGridPartIterator< GridPart, codim > Iterator;

      static const std::string itName = std::string( rgName ) + "Iterator";
      pybind11::class_< Iterator > itCls( scope, itName.c_str() );
      itCls.def( "__iter__", [] ( Iterator &it ) -> Iterator & { return it; } );
      itCls.def( "__next__", &Iterator::next );

      pybind11::class_< Range > rgCls( scope, rgName );
      rgCls.def( "__iter__", [] ( const Range &range ) { return Iterator( range ); } );
    };




    // makePyGlobalGridFunction
    // ------------------------

    template< class GridPart, int dimRange >
    auto makePyGlobalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
    {
      typedef typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate Coordinate;
      return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Coordinate &x ) {
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate.call( x ) );
          return v.template cast< FieldVector< double, dimRange > >();
        } );
    }



    // makePyLocalGridFunction
    // -----------------------

    template< class GridPart, int dimRange >
    auto makePyLocalGridFunction ( const GridPart &gridPart, std::string name, pybind11::function evaluate, std::integral_constant< int, dimRange > )
    {
      typedef typename GridPart::template Codim< 0 >::EntityType Entity;
      typedef typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate Coordinate;
      return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Entity &entity, const Coordinate &x ) {
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate.call( entity, x ) );
          return v.template cast< FieldVector< double, dimRange > >();
        } );
    }



    // registerPyGlobalGridFunction
    // ----------------------------

    template< class GridPart, int dimRange >
    void registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
    {
      typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
      static const std::string clsName = name + std::to_string( dimRange );
      registerGridFunction< GridFunction >( scope, clsName.c_str() );
    };

    template< class GridPart, int... dimRange >
    void registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( (registerPyGlobalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() ), dimRange)... );
    };



    // registerPyLocalGridFunction
    // ---------------------------

    template< class GridPart, int dimRange >
    void registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
    {
      typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
      static const std::string clsName = name + std::to_string( dimRange );
      registerGridFunction< GridFunction >( scope, clsName.c_str() );
    };

    template< class GridPart, int... dimRange >
    void registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( (registerPyLocalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() ), dimRange)... );
    };



    // arrayFromIntegerSequence
    // ------------------------

    template< class T, class Make, class I, I... i, class... Args >
    std::array< T, sizeof...( i ) > arrayFromIntegerSequence ( Make make, std::integer_sequence< I, i... >, const Args &... args )
    {
      return {{ T( make( std::integral_constant< I, i >(), args... ) )... }};
    }



    // registerGrid
    // ------------

    template< class GridPart >
    void registerGrid ( pybind11::module module )
    {
      const int dim = GridPart::dimension;

      registerHierarchicalGrid< typename GridPart::GridType >( module );

      typedef LeafGrid< GridPart > G;
      pybind11::class_< G > grid( module, "LeafGrid" );
      grid.def( pybind11::init< std::string >() );

      registerPyGridPartRange< GridPart, 0 >( grid, "Elements" );
      grid.def_property_readonly( "elements", [] ( pybind11::object g ) {
          return PyGridPartRange< GridPart, 0 >( *g.template cast< const G & >().gridPart(), g );
        } );

      registerPyGridPartRange< GridPart, dim >( grid, "Vertices" );
      grid.def_property_readonly( "vertices", [] ( pybind11::object g ) {
          return PyGridPartRange< GridPart, dim >( *g.template cast< const G & >().gridPart(), g );
        } );

      grid.def( "__repr__", [] ( const G &grid ) -> std::string {
          return "LeafGrid with " + std::to_string( grid.size( 0 ) ) + " elements";
        } );

      grid.def_property_readonly( "hierarchicalGrid", &G::hierarchicalGrid );

      grid.def( "size", &G::size );

      const int maxDimRange = 10;
      registerPyGlobalGridFunction< GridPart >( grid, "GlobalGridFunction", std::make_integer_sequence< int, maxDimRange+1 >() );
      registerPyLocalGridFunction< GridPart >( grid, "LocalGridFunction", std::make_integer_sequence< int, maxDimRange+1 >() );

      typedef std::function< pybind11::object( const GridPart &, std::string, pybind11::function, pybind11::object ) > DispatchGridFunction;
      auto makeDispatchGlobalGridFunction = [] ( auto dimRange ) {
          return [ dimRange ] ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent ) {
              return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), std::move( evaluate ), dimRange ), pybind11::return_value_policy::move, parent );
            };
        };
      auto dispatchGlobalGridFunction = arrayFromIntegerSequence< DispatchGridFunction >( makeDispatchGlobalGridFunction, std::make_integer_sequence< int, maxDimRange+1 >() );

      grid.def( "globalGridFunction", [ dispatchGlobalGridFunction ] ( pybind11::object g, std::string name, pybind11::function evaluate ) {
          const GridPart &gridPart = *g.cast< const G & >().gridPart();
          typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate x( 0 );
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate.call( x ) );
          const int dimRange = len( v );
          if( (dimRange > maxDimRange) )
            DUNE_THROW( NotImplemented, "globalGridFunction not implemented for dimRange = " + std::to_string( dimRange ) );
          return dispatchGlobalGridFunction[ dimRange ]( gridPart, std::move( name ), std::move( evaluate ), std::move( g ) );
        } );

      auto makeDispatchLocalGridFunction = [] ( auto dimRange ) {
          return [ dimRange ] ( const GridPart &gridPart, std::string name, pybind11::function evaluate, pybind11::object parent ) {
              return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), std::move( evaluate ), dimRange ), pybind11::return_value_policy::move, parent );
            };
        };
      auto dispatchLocalGridFunction = arrayFromIntegerSequence< DispatchGridFunction >( makeDispatchLocalGridFunction, std::make_integer_sequence< int, maxDimRange+1 >() );

      grid.def( "localGridFunction", [ dispatchLocalGridFunction ] ( pybind11::object g, std::string name, pybind11::function evaluate ) {
          const GridPart &gridPart = *g.cast< const G & >().gridPart();
          int dimRange = -1;
          if( gridPart.template begin< 0 >() != gridPart.template end< 0 >() )
          {
            typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate x( 0 );
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate.call( *gridPart.template begin< 0 >(), x ) );
            dimRange = len( v );
          }
          dimRange = gridPart.comm().max( dimRange );
          if( dimRange < 0 )
            DUNE_THROW( InvalidStateException, "Cannot create local grid function on empty grid" );
          if( dimRange > maxDimRange )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimRange ) );
          return dispatchLocalGridFunction[ dimRange ]( gridPart, std::move( name ), std::move( evaluate ), std::move( g ) );
        } );

      registerGridFunctionInterface< GridPart >( module, Std::make_integer_sequence< int, 10 >() );
//      registerGridFunctionExpression< GridPart >( Std::make_integer_sequence< int, 10 >() );
//      registerLocalGridFunctionExpression< GridPart >( Std::make_integer_sequence< int, 10 >() );

      registerVTKWriter< VTKWriter< typename GridPart::GridViewType > >( module );
      grid.def( "vtkWriter", [] ( const G &grid ) {
          return new VTKWriter< typename GridPart::GridViewType >( static_cast< typename GridPart::GridViewType >( *grid.gridPart() ) );
        }, pybind11::keep_alive< 0, 1 >() );

      typedef VTKOutput< GridPart > VTK;

      void (VTK::*add1)( std::shared_ptr< GridFunction<GridPart,1> > ) = &VTK::template add<1>;
      void (VTK::*add2)( std::shared_ptr< GridFunction<GridPart,2> > ) = &VTK::template add<2>;
      void (VTK::*add3)( std::shared_ptr< GridFunction<GridPart,3> > ) = &VTK::template add<3>;
      void (VTK::*add4)( std::shared_ptr< GridFunction<GridPart,4> > ) = &VTK::template add<4>;
      void (VTK::*add5)( std::shared_ptr< GridFunction<GridPart,5> > ) = &VTK::template add<5>;
      void (VTK::*add6)( std::shared_ptr< GridFunction<GridPart,6> > ) = &VTK::template add<6>;
      void (VTK::*add7)( std::shared_ptr< GridFunction<GridPart,7> > ) = &VTK::template add<7>;
      void (VTK::*add8)( std::shared_ptr< GridFunction<GridPart,8> > ) = &VTK::template add<8>;
      void (VTK::*add9)( std::shared_ptr< GridFunction<GridPart,9> > ) = &VTK::template add<9>;
      pybind11::class_< VTK > vtk( module, "VTKOutput" );
      vtk.def( pybind11::init< G >() );
      vtk.def( "write", &VTK::write );
      vtk.def( "add", add1 );
      vtk.def( "add", add2 );
      vtk.def( "add", add3 );
      vtk.def( "add", add4 );
      vtk.def( "add", add5 );
      vtk.def( "add", add6 );
      vtk.def( "add", add7 );
      vtk.def( "add", add8 );
      vtk.def( "add", add9 );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PYGRID_HH
