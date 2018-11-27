#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/fem/misc/domainintegral.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/grid/vtk.hh>

#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/gridfunctionview.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/numpy.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/pybind11/pybind11.hh>

namespace Dune
{

  namespace FemPy
  {

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction, class... options >
    inline static void registerLocalFunction ( pybind11::handle scope, pybind11::class_< LocalFunction, options... > cls)
    {
      typedef typename LocalFunction::EntityType::Geometry::LocalCoordinate LocalCoordinate;
      Dune::Python::registerFieldVector<double,LocalFunction::RangeType::dimension>(scope);
      Dune::Python::registerFieldMatrix<double, LocalFunction::RangeType::dimension, LocalFunction::EntityType::Geometry::coorddimension>( scope );
      cls.def_property_readonly( "dimRange", [] ( LocalFunction & ) -> int { return LocalFunction::RangeType::dimension; } );
      cls.def_property_readonly( "order", [] ( const LocalFunction & lf ) -> int { return lf.order(); } );
      cls.def( "evaluate", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          lf.evaluate( x, value );
          return value;
        } );
      cls.def( "jacobian", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::JacobianRangeType jacobian;
          lf.jacobian( x, jacobian );
          return jacobian;
        } );
      cls.def( "hessian", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::HessianRangeType hessian;
          lf.hessian( x, hessian );
          return hessian;
        } );
    }



    namespace detail
    {

      // makePyLocalGridFunction
      // -----------------------

      template< class GridPart, int dimRange >
      inline static auto makePyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::EntityType Entity;
        typedef typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Entity &entity, const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( entity, x ) );
            try { return v.template cast< FieldVector< double, dimRange > >(); }
            catch(std::exception& e) { std::cout << e.what() << " in LocalGridFunction::evaluate" << std::endl; throw pybind11::cast_error("error converting return value in localGridFunction"); }
            return FieldVector<double,dimRange>(0);
          }, order );
      }



      template< class GridPart, int dimRange >
      auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > );




      // registerGridFunctionName
      // ------------------------

      template< class GridFunction, class... options >
      inline static auto registerGridFunctionName ( pybind11::class_< GridFunction, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_same< decltype( std::declval< GridFunction & >().name() ), std::string & >::value >
      {
        cls.def_property( "name", [] ( GridFunction &self ) -> std::string { return self.name(); }, [] ( GridFunction &self, std::string name ) { self.name() = name; } );
      }

      template< class GridFunction, class... options >
      inline static void registerGridFunctionName ( pybind11::class_< GridFunction, options... > cls, PriorityTag< 0 > )
      {
        cls.def_property_readonly( "name", [] ( GridFunction &self ) -> std::string { return self.name(); } );
      }

      template< class GridFunction, class... options >
      inline static void registerGridFunctionName ( pybind11::class_< GridFunction, options... > cls )
      {
        registerGridFunctionName( cls, PriorityTag< 42 >() );
      }



      // registerGridFunction
      // --------------------

      template< class GridFunction, class... options >
      inline static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
      {
        using pybind11::operator""_a;

        // typedef typename GridFunction::LocalFunctionType LocalFunction;
        typedef typename Dune::Fem::ConstLocalFunction<GridFunction> LocalFunction;
        typedef typename LocalFunction::EntityType Entity;
        typedef typename GridFunction::GridPartType GridPartType;
        typedef typename GridPartType::GridViewType GridView;

        auto lfClass = Python::insertClass<LocalFunction>(cls, "LocalFunction",
            Python::GenerateTypeName("TODO-LGF"));
            // Python::GenerateTypeName(cls,"LocalFunctionType")).first;
        assert( lfClass.second );
        registerLocalFunction< LocalFunction >( cls, lfClass.first );

        /*
        cls.def( "__repr__", [] ( GridFunction &self ) -> std::string {
            return "GridFunction< " + std::to_string( GridFunction::RangeType::dimension ) + " >(name = " + self.name() + ")";
          } );
        */

        cls.def_property_readonly( "dimRange", [] ( GridFunction & ) -> int { return GridFunction::RangeType::dimension; } );
        cls.def_property_readonly( "order", [] ( GridFunction &self ) -> int { return self.space().order(); } );
        cls.def_property_readonly( "grid", [] ( GridFunction &self ) -> GridView { return static_cast< GridView >( self.gridPart() ); } );
        cls.def_property_readonly( "space", &GridFunction::space );

        registerGridFunctionName( cls );

        cls.def( "localFunction", [] ( const GridFunction &self, const Entity &entity ) { // -> LocalFunction {
            auto ret = std::make_unique<LocalFunction>(self);
            ret->bind(entity);
            return ret; // using ConstLocalFunction now instead of self.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >(),
             pybind11::return_value_policy::take_ownership );

        typedef decltype( makePyLocalGridFunction( std::declval< GridPartType >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, 1 >() ) ) ScalarLocalGridFunction;
        if (!pybind11::already_registered<ScalarLocalGridFunction>())
          registerPyLocalGridFunction<GridPartType>( scope, "LocalGridFunction", std::integral_constant< int, 1 >() );

        cls.def( "__getitem__", [] ( const GridFunction &self, std::size_t c ) {
            return makePyLocalGridFunction( self.gridPart(), self.name() + "_" + std::to_string(c), self.space().order(),
                pybind11::cpp_function( [ self, c ] ( const Entity &e, const typename Entity::Geometry::LocalCoordinate &x ) {
                    return Fem::ConstLocalFunction<GridFunction>(e,self).evaluate(x);
                  } ), std::integral_constant< int, 1 >() );
          }, pybind11::keep_alive< 0, 1 >() );

        cls.def( "addToVTKWriter", &Dune::Python::addToVTKWriter< GridFunction >, pybind11::keep_alive< 3, 1 >(), "name"_a, "writer"_a, "dataType"_a );

        cls.def( "cellData", [] ( const GridFunction &self, int level ) { return cellData( self, refinementLevels( level ) ); }, "level"_a = 0 );
        cls.def( "pointData", [] ( const GridFunction &self, int level ) { return pointData( self, refinementLevels( level ) ); }, "level"_a = 0 );

        cls.def( "integrate", [] ( const GridFunction &self ) { return Dune::Fem::Integral< GridPartType >( self.gridPart(), self.space().order() ).norm( self ); } );


#if 0
        cls.def_property_readonly( "as_ufl", [] ( GridFunction &self ) -> pybind11::handle {
              pybind11::tuple args( 1 );
              args[ 0 ] = self;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            } );
#elif 0
        cls.def( "as_ufl", [] ( GridFunction &self ) -> pybind11::handle {
              pybind11::tuple args( 1 );
              args[ 0 ] = self;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            }, pybind11::keep_alive<0,1>() );
#else
        cls.def( "as_ufl", [] ( pybind11::object &self ) -> pybind11::handle {
              pybind11::tuple args( 1 );
              args[ 0 ] = self;
              return PyObject_Call( Dune::FemPy::getGridFunctionWrapper().ptr(), args.ptr(), nullptr );
            },  pybind11::keep_alive< 0, 1 >() );
#endif
      }



      // clsVirtualizedGridFunction
      // --------------------------

      template< class GridPart, class Value >
      inline pybind11::class_< VirtualizedGridFunction< GridPart, Value > >
      clsVirtualizedGridFunction ( pybind11::handle scope )
      {
        typedef VirtualizedGridFunction< GridPart, Value > GridFunction;
        std::string adapt = Python::GenerateTypeName("Dune::FemPy::GridPartAdapter",
            Dune::MetaType<typename GridPart::GridViewType>()).name();
        auto vgfClass = Python::insertClass<GridFunction>(scope,"VirtualizedGridFunction"+std::to_string(Value::dimension),
            Python::GenerateTypeName("VirtualizedGridFunction",adapt,MetaType<Value>()),
            Python::IncludeFiles{"dune/fempy/function/virtualizedgridfunction.hh"});

        if( vgfClass.second )
          detail::registerGridFunction( scope, vgfClass.first );

        typedef typename GridFunction::DiscreteFunctionSpaceType SpaceType;
        if( !pybind11::already_registered< SpaceType >() )
        {
          auto spcCls = Python::insertClass<SpaceType>(vgfClass.first, "Space",
              Python::GenerateTypeName("TODO-SPACE"));
          if (spcCls.second)
            detail::registerFunctionSpace( vgfClass.first, spcCls.first );
        }

        return vgfClass.first;
      }

    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction, class... options >
    inline static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Value;

      detail::registerGridFunction( scope, cls );
      detail::clsVirtualizedGridFunction< GridPart, Value >( scope ).def( pybind11::init< GridFunction >() );
      pybind11::implicitly_convertible< GridFunction, VirtualizedGridFunction< GridPart, Value > >();
    }

    // registerVirtualizedGridFunction
    // -------------------------------

    template< class GridPart, int... dimRange >
    inline static void registerVirtualizedGridFunction ( pybind11::handle scope, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::clsVirtualizedGridFunction< GridPart, FieldVector< double, dimRange > >( scope )... );
    };



    namespace detail
    {

      // makePyGlobalGridFunction
      // ------------------------

      template< class GridPart, int dimRange >
      inline static auto makePyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, std::integral_constant< int, dimRange > )
      {
        typedef typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate Coordinate;
        return simpleGridFunction( std::move( name ), gridPart, [ evaluate ] ( const Coordinate &x ) {
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( x ) );
            try { return v.template cast< FieldVector< double, dimRange > >(); }
            catch(std::exception& e) { std::cout << e.what() << " in GlobalGridFunction::evaluate" << std::endl; throw pybind11::cast_error("error converting return value in localGridFunction"); }
            return FieldVector<double,dimRange>(0);
          }, order );
      }



      // registerPyGlobalGridFunction
      // ----------------------------

      template< class GridPart, int dimRange >
      inline static auto registerPyGlobalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyGlobalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        std::string adapt = Python::GenerateTypeName("Dune::FemPy::GridPartAdapter",
            Dune::MetaType<typename GridPart::GridViewType>()).name();
        auto gfClass = Python::insertClass<GridFunction>(scope, name+std::to_string(dimRange),
            Python::GenerateTypeName("Dune::FemPy::VirtualizeGridFunction",
              adapt,
              "Dune::FieldVector<double,"+std::to_string(dimRange)+">"),
              Python::IncludeFiles{"dune/fempy/function/virtualizedgridfunction.hh","dune/fempy/grid/gridpartadapter.hh"}
              );
        if (gfClass.second)
          FemPy::registerGridFunction< GridFunction >( scope, gfClass.first );
        typedef typename GridFunction::DiscreteFunctionSpaceType SpaceType;
        if( !pybind11::already_registered< SpaceType >() )
        {
          auto spcCls = Python::insertClass<SpaceType>(gfClass.first, "Space",
              Python::GenerateTypeName("TODO-SPACE"));
          if (spcCls.second)
            detail::registerFunctionSpace( gfClass.first, spcCls.first );
        }

        return gfClass.first;
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      inline static pybind11::object pyGlobalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyGlobalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defGlobalGridFunction
    // ---------------------

    template< class GridPart, int... dimRange >
    inline static auto defGlobalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyGlobalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, int, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyGlobalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gv, std::string name, int order, pybind11::function evaluate ) {
          typedef typename GridPart::GridViewType GridView;
          const auto &gp = gridPart<GridView>( gv );
          // typename GridPart::template Codim< 0 >::GeometryType::GlobalCoordinate x( 0 );
          auto x = gp.template begin<0>()->geometry().center();
          pybind11::gil_scoped_acquire acq;
          pybind11::object v( evaluate( x ) );
          const std::size_t dimR = len( v );
          if( dimR >= dispatch.size() )
            DUNE_THROW( NotImplemented, "globalGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ dimR ]( gp, std::move( name ), order, std::move( evaluate ), std::move(gv) );
        };
    }



    namespace detail
    {

      // registerPyLocalGridFunction
      // ---------------------------

      template< class GridPart, int dimRange >
      inline auto registerPyLocalGridFunction ( pybind11::handle scope, const std::string &name, std::integral_constant< int, dimRange > )
      {
        typedef decltype( makePyLocalGridFunction( std::declval< GridPart >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, dimRange >() ) ) GridFunction;
        const std::string clsName = name + std::to_string( dimRange );
        std::string adapt = Python::GenerateTypeName("Dune::FemPy::GridPartAdapter",
            Dune::MetaType<typename GridPart::GridViewType>()).name();
        auto gfClass = Python::insertClass<GridFunction>(scope,name+std::to_string(dimRange),
            Python::GenerateTypeName("Dune::FemPy::VirtualizeGridFunction",
              adapt,
              "Dune::FieldVector<double,"+std::to_string(dimRange)+">"),
              Python::IncludeFiles{"dune/fempy/function/virtualizedgridfunction.hh","dune/fempy/grid/gridpartadapter.hh"}
              );

        if (gfClass.second)
          FemPy::registerGridFunction( scope, gfClass.first);

        typedef typename GridFunction::DiscreteFunctionSpaceType SpaceType;
        if( !pybind11::already_registered< SpaceType >() )
        {
          auto spcCls = Python::insertClass<SpaceType>(gfClass.first, "Space",
              Python::GenerateTypeName("TODO-SPACE"));
          if (spcCls.second)
            detail::registerFunctionSpace( gfClass.first, spcCls.first );
        }

        return gfClass.first;
      };



      // pyGlobalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      inline static pybind11::object pyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail



    // defLocalGridFunction
    // --------------------

    template< class GridPart, int... dimRange >
    inline static auto defLocalGridFunction ( pybind11::handle scope, std::string name, std::integer_sequence< int, dimRange... > )
    {
      std::ignore = std::make_tuple( detail::registerPyLocalGridFunction< GridPart >( scope, name, std::integral_constant< int, dimRange >() )... );

      typedef std::function< pybind11::object( const GridPart &, std::string, int, pybind11::function, pybind11::object ) > Dispatch;
      std::array< Dispatch, sizeof...( dimRange ) > dispatch = {{ Dispatch( detail::pyLocalGridFunction< GridPart, dimRange > )... }};

      return [ dispatch ] ( pybind11::object gv, std::string name, int order, pybind11::function evaluate ) {
          typedef typename GridPart::GridViewType GridView;
          const auto &gp = gridPart<GridView>( gv );
          int dimR = -1;
          if( gp.template begin< 0 >() != gp.template end< 0 >() )
          {
            typename GridPart::template Codim< 0 >::GeometryType::LocalCoordinate x( 0 );
            pybind11::gil_scoped_acquire acq;
            pybind11::object v( evaluate( *gp.template begin< 0 >(), x ) );
            dimR = len( v );
          }
          dimR = gp.comm().max( dimR );
          if( dimR < 0 )
            DUNE_THROW( InvalidStateException, "Cannot create local grid function on empty grid" );
          if( static_cast< std::size_t >( dimR ) >= dispatch.size() )
            DUNE_THROW( NotImplemented, "localGridFunction not implemented for dimRange = " + std::to_string( dimR ) );
          return dispatch[ static_cast< std::size_t >( dimR ) ]( gp, std::move( name ), order, std::move( evaluate ), std::move( gv ) );
        };
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
