#ifndef DUNE_FEMPY_PYGRID_HH
#define DUNE_FEMPY_PYGRID_HH

#include <list>
#include <memory>
#include <tuple>

#include <dune/common/std/utility.hh>

#include <dune/fempy/grid.hh>
#include <dune/fempy/pybind11/functional.h>
#include <dune/fempy/pybind11/pybind11.h>
#include <dune/fempy/pybind11/stl.h>
#include <dune/fempy/pygridfunction.hh>

namespace Dune
{

  namespace FemPy
  {

    // registerIteratorRange
    // ---------------------

    template< class IteratorRange >
    void registerIteratorRange ( pybind11::handle scope, const char *name )
    {
      pybind11::class_< IteratorRange > cls( scope, name );
      cls.def( "__iter__", [] ( const IteratorRange &rg ) { return pybind11::make_iterator( rg.begin(), rg.end() ); }, pybind11::keep_alive< 0, 1 >() );
    }



    // registerGridGeometry
    // --------------------

    template< class Geometry >
    void registerGridGeometry ( pybind11::module module )
    {
      typedef typename Geometry::LocalCoordinate LocalCoordinate;
      typedef typename Geometry::ctype ctype;

      typedef FemPy::CornerIterator< Geometry > CornerIterator;

      static const std::string clsName = "Geometry" + std::to_string( Geometry::mydimension );
      pybind11::class_< Geometry > cls( module, clsName.c_str() );

      registerIteratorRange< IteratorRange< CornerIterator > >( cls, "Corners" );
      cls.def_property_readonly( "corners", [] ( const Geometry &geo ) {
          return IteratorRange< CornerIterator >( CornerIterator( geo, 0 ), CornerIterator( geo ) );
        }, pybind11::keep_alive< 0, 1 >() );

      cls.def_property_readonly( "center", &Geometry::center );
      cls.def_property_readonly( "volume", &Geometry::volume );

      cls.def_property_readonly( "affine", &Geometry::affine );

      cls.def( "global", &Geometry::global );
      cls.def( "global", [] ( const Geometry &geo, const pybind11::list &x ) {
          LocalCoordinate y( 0 );
          const std::size_t size = x.size();
          for( std::size_t i = 0; i < size; ++i )
            y[ i ] = x[ i ].template cast< ctype >();
          return geo.global( y );
        } );

      cls.def( "integrationElement", &Geometry::integrationElement );
      cls.def( "integrationElement", [] ( const Geometry &geo, const pybind11::list &x ) {
          LocalCoordinate y( 0 );
          const std::size_t size = x.size();
          for( std::size_t i = 0; i < size; ++i )
            y[ i ] = x[ i ].template cast< ctype >();
          return geo.integrationElement( y );
        } );
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



    // registerGrid
    // ------------

    template< class GridPart >
    void registerGrid ( pybind11::module module )
    {
      registerHierarchicalGrid< typename GridPart::GridType >( module );

      typedef LeafGrid< GridPart > G;
      pybind11::class_< G > grid( module, "LeafGrid" );
      grid.def( pybind11::init< std::string >() );

      registerIteratorRange< decltype( elements( std::declval< GridPart >(), Partitions::interiorBorder ) ) >( grid, "Elements" );
      grid.def_property_readonly( "elements", [] ( const G &grid ) { return elements( *grid.gridPart(), Partitions::interiorBorder ); }, pybind11::keep_alive< 0, 1 >() );

      registerIteratorRange< decltype( vertices( std::declval< GridPart >(), Partitions::interiorBorder ) ) >( grid, "Vertices" );
      grid.def_property_readonly( "vertices", [] ( const G &grid ) { return vertices( *grid.gridPart(), Partitions::interiorBorder ); }, pybind11::keep_alive< 0, 1 >() );

      grid.def( "__repr__", [] ( const G &grid ) -> std::string {
          return "LeafGrid with " + std::to_string( grid.size( 0 ) ) + " elements";
        } );

      grid.def_property_readonly( "hierarchicalGrid", &G::hierarchicalGrid );

      grid.def( "size", &G::size );

      registerGridFunctionInterface< GridPart >( module, Std::make_integer_sequence< int, 10 >() );
//      registerGridFunctionExpression< GridPart >( Std::make_integer_sequence< int, 10 >() );
//      registerLocalGridFunctionExpression< GridPart >( Std::make_integer_sequence< int, 10 >() );

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
