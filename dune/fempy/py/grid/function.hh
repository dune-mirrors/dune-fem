#ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
#define DUNE_FEMPY_PY_GRID_FUNCTION_HH

#include <functional>
#include <string>
#include <tuple>
#include <utility>

#include <dune/common/typeutilities.hh>
#include <dune/common/visibility.hh>

#include <dune/fem/misc/domainintegral.hh>
#include <dune/fem/misc/gridfunctionview.hh>

#include <dune/python/common/typeregistry.hh>
#include <dune/python/grid/vtk.hh>
#include <dune/python/grid/localview.hh>

#include <dune/fempy/function/simplegridfunction.hh>
#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/grid/gridpart.hh>
#include <dune/fempy/py/grid/numpy.hh>
#include <dune/fempy/py/space.hh>
#include <dune/fempy/py/grid/gridpart.hh>
// #include <dune/fempy/py/function/grid.hh>
#include <dune/fempy/pybind11/pybind11.hh>

#if HAVE_DUNE_VTK
#include <dune/vtk/function.hh>
#endif

namespace Dune
{

  namespace FemPy
  {

    // registerLocalFunction
    // ---------------------

    template< class LocalFunction, class... options >
    inline static void registerLocalFunction ( pybind11::handle scope, pybind11::class_< LocalFunction, options... > cls)
    {
      typedef typename LocalFunction::EntityType Element;
      typedef typename Element::Geometry::LocalCoordinate LocalCoordinate;
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
      cls.def( "__call__", [] ( const LocalFunction &lf, const LocalCoordinate &x ) {
          typename LocalFunction::RangeType value;
          lf.evaluate( x, value );
          bool scalar = pybind11::cast(lf).attr("scalar").template cast<bool>();
          return scalar? pybind11::cast(value[0]) : pybind11::cast(value);
        } );
      Dune::Python::registerLocalView< Element >( cls );
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

        typedef typename Dune::Fem::ConstLocalFunction<GridFunction> LocalFunction;
        typedef typename LocalFunction::EntityType Entity;
        typedef typename GridFunction::GridPartType GridPartType;
        typedef typename GridPartType::GridViewType GridView;

        auto lfClass = Python::insertClass<LocalFunction>(cls, "LocalFunction",
            Python::GenerateTypeName("TODO-LGF"), pybind11::dynamic_attr());
        assert( lfClass.second );
        registerLocalFunction< LocalFunction >( cls, lfClass.first );

        cls.def_property_readonly( "order", [] ( GridFunction &self ) -> int { return self.space().order(); } );
        cls.def_property_readonly( "grid", [] ( GridFunction &self ) -> const GridView&
          {
            PyErr_WarnEx(PyExc_DeprecationWarning, "attribute 'grid' is deprecated, use 'gridView' instead.", 2);
            return self.gridPart();
          } );
        cls.def_property_readonly( "gridView", [] ( GridFunction &self ) -> const GridView& { return self.gridPart(); } );

        registerGridFunctionName( cls );

        cls.def( "localFunction", [] ( const GridFunction &self ) { // -> LocalFunction {
            auto ret = std::make_unique<LocalFunction>(self);
            return ret;
          }, pybind11::keep_alive< 0, 1 >(),
             pybind11::return_value_policy::take_ownership );
        cls.def( "localFunction", [] ( const GridFunction &self, const Entity &entity ) { // -> LocalFunction {
            auto ret = std::make_unique<LocalFunction>(self);
            ret->bind(entity);
            return ret; // using ConstLocalFunction now instead of self.localFunction( entity );
          }, pybind11::keep_alive< 0, 1 >(), pybind11::keep_alive< 0, 2 >(),
             pybind11::return_value_policy::take_ownership );

        typedef decltype( makePyLocalGridFunction( std::declval< GridPartType >(), std::declval< std::string >(), std::declval< int >(), std::declval< pybind11::function >(), std::integral_constant< int, 1 >() ) ) ScalarLocalGridFunction;
        if (!pybind11::already_registered<ScalarLocalGridFunction>())
          registerPyLocalGridFunction<GridPartType>( cls, "LocalGridFunction", std::integral_constant< int, 1 >() );

        cls.def( "__getitem__", [] ( const GridFunction &self, std::size_t c ) {
            const auto pself = &self; // make sure reference is correctly capture in lambda and doesn't go out of scope
            return makePyLocalGridFunction( self.gridPart(), self.name() + "_" + std::to_string(c), self.space().order(),
                pybind11::cpp_function( [ pself, c ] ( const Entity &e, const typename Entity::Geometry::LocalCoordinate &x ) {
                    return Fem::ConstLocalFunction<GridFunction>(e,*pself).evaluate(x)[c];
                  } ), std::integral_constant< int, 1 >() );
          }, pybind11::keep_alive< 0, 1 >() );

        cls.def( "integrate", [] ( pybind11::handle self )
            { const GridFunction &gf = self.template cast<GridFunction>();
              auto value = Dune::Fem::Integral< GridPartType >( gf.gridPart(), gf.space().order() ).norm( gf );
              return self.attr("scalar").template cast<bool>() ?
                 pybind11::cast(value[0]) : pybind11::cast(value);
            } );
        cls.def( "as_ufl", [] ( pybind11::object self ) {
              return Dune::FemPy::getGridFunctionWrapper()(self);
            } );

        Dune::Python::detail::registerGridFunction(scope,cls);
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
        {
          Dune::FemPy::detail::registerGridFunction( scope, vgfClass.first );
        }

        return vgfClass.first;
      }

    } // namespace detail



    // registerGridFunction
    // --------------------

    template< class GridFunction >
    inline static void registerGridFunction ( pybind11::handle scope )
    {
      typedef typename GridFunction::GridView GridView;
      typedef typename GridFunction::Value Value;
      detail::clsVirtualizedGridFunction< GridView, Value >( scope ).
             def( pybind11::init< GridFunction >() );
      pybind11::implicitly_convertible< GridFunction, VirtualizedGridFunction< GridView, Value > >();
#if HAVE_DUNE_VTK
      using GridView = typename GridFunction::GridView;
      using VtkGF = Dune::Vtk::Function<GridView>;
      // register the Function class if not already available
      auto vgfClass = Python::insertClass<VtkGF>(scope,"VtkFunction",
          Python::GenerateTypeName("Dune::Vtk::Function",MetaType<GridView>()),
          Python::IncludeFiles{"dune/vtk/function.hh"});
      assert( !vgfClass.second );
      vgfClass.first.def( pybind11::init( [] ( GridFunction &gf ) {
          return new VtkGF( localFunction(gf), gf.name() );
        } ) );
      pybind11::implicitly_convertible<GridFunction,VtkGF>();
#endif
    }

    template< class GridFunction, class... options >
    inline static void registerGridFunction ( pybind11::handle scope, pybind11::class_< GridFunction, options... > cls )
    {
      typedef typename GridFunction::GridPartType GridPart;
      typedef typename GridFunction::RangeType Value;

      // testing with uflFunction overal 36s:  18s base time before getting here
      // 13s for first line, 4s second line
      Dune::FemPy::detail::registerGridFunction( scope, cls );
      detail::clsVirtualizedGridFunction< GridPart, Value >( scope ).def( pybind11::init< GridFunction >() );
      pybind11::implicitly_convertible< GridFunction, VirtualizedGridFunction< GridPart, Value > >();

#if HAVE_DUNE_VTK
      using GridView = typename GridFunction::GridView;
      using VtkGF = Dune::Vtk::Function<GridView>;
      // register the Function class if not already available
      auto vgfClass = Python::insertClass<VtkGF>(scope,"VtkFunction",
          Python::GenerateTypeName("Dune::Vtk::Function",MetaType<GridView>()),
          Python::IncludeFiles{"dune/vtk/function.hh"});
      assert( !vgfClass.second );
      vgfClass.first.def( pybind11::init( [] ( GridFunction &gf ) {
          return new VtkGF( localFunction(gf), gf.name() );
        } ) );
      pybind11::implicitly_convertible<GridFunction,VtkGF>();
#endif
    }


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
              Python::IncludeFiles{"dune/fempy/function/virtualizedgridfunction.hh","dune/fem/gridpart/common/gridpartadapter.hh"}
              );
        if (gfClass.second)
        {
          Dune::FemPy::registerGridFunction< GridFunction >( scope, gfClass.first );
          gfClass.first.def_property_readonly( "scalar", [] ( GridFunction &self ) { return false; } );
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
              Python::IncludeFiles{"dune/fempy/function/virtualizedgridfunction.hh","dune/fem/gridpart/common/gridpartadapter.hh"}
              );

        if (gfClass.second)
        {
          Dune::FemPy::registerGridFunction( scope, gfClass.first);
          gfClass.first.def_property_readonly( "scalar", [] ( GridFunction &self ) { return false; } );
        }
        return gfClass.first;
      };



      // pyLocalGridFunction
      // --------------------

      template< class GridPart, int dimRange >
      inline static pybind11::object pyLocalGridFunction ( const GridPart &gridPart, std::string name, int order, pybind11::function evaluate, pybind11::object parent )
      {
        return pybind11::cast( makePyLocalGridFunction( gridPart, std::move( name ), order, std::move( evaluate ), std::integral_constant< int, dimRange >() ), pybind11::return_value_policy::move, parent );
      }

    } // namespace detail

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_GRID_FUNCTION_HH
