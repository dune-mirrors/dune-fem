#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fempy/pybind11/pybind11.hh>

#include <dune/common/typeutilities.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/bcrsmatrix.hh>
#include <dune/corepy/istl/bcrsmatrix.hh>
#endif // #if HAVE_DUNE_ISTL

#include <dune/fem/misc/l2norm.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/py/space.hh>


namespace Dune
{

  namespace FemPy
  {

    // registerScheme
    // --------------

    namespace detail
    {

#if HAVE_DUNE_ISTL
      template< class B, class A >
      inline static const BCRSMatrix< B, A > &getBCRSMatrix ( const BCRSMatrix< B, A > &matrix ) noexcept
      {
        return matrix;
      }
#endif // #if HAVE_DUNE_ISTL



      // registerSchemeConstructor
      // -------------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> std::enable_if_t< std::is_constructible< Scheme, const typename Scheme::DiscreteFunctionSpaceType &, const typename Scheme::ModelType & >::value >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::ModelType ModelType;

        using pybind11::operator""_a;

        cls.def( pybind11::init( [] ( Space &space, const ModelType &model ) {
            return new Scheme( space, model );
          } ), "space"_a, "model"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
        cls.def( pybind11::init( [] ( Space &space, const ModelType &model, const pybind11::dict &parameters ) {
            return new Scheme( space, model, pyParameter( parameters, std::make_shared< std::string >() ) );
          } ), "space"_a, "model"_a, "parameters"_a, pybind11::keep_alive< 1, 2 >(), pybind11::keep_alive< 1, 3 >() );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeConstructor ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeConstructor( cls, PriorityTag< 42 >() );
      }



      // registerSchemeAssemble
      // ----------------------

      // register assemble method if data method is available (and return value is registered)
#if HAVE_DUNE_ISTL
      template< class Scheme, class... options >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 2 > )
        -> void_t< decltype( getBCRSMatrix( std::declval< const typename Scheme::LinearOperatorType & >().matrix() ) ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        typedef std::decay_t< decltype( getBCRSMatrix( std::declval< const typename Scheme::LinearOperatorType & >().matrix() ) ) > BCRSMatrix;
        if( !pybind11::already_registered< BCRSMatrix >() )
          CorePy::registerBCRSMatrix< BCRSMatrix >( cls );

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( Scheme &self, const DiscreteFunction &ubar ) {
            return getBCRSMatrix( self.assemble( ubar ).matrix() );
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
        cls.def( "assemble", [] ( Scheme &self, const VirtualizedGridFunction< GridPart, RangeType > &ubar ) {
            return getBCRSMatrix( self.assemble( ubar ).matrix() );
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
      }
#endif // #if HAVE_DUNE_ISTL

      template< class Scheme, class... options >
      inline static auto registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< const typename Scheme::LinearOperatorType & >().matrix().data() ) >
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename DiscreteFunction::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        using pybind11::operator""_a;

        cls.def( "assemble", [] ( Scheme &self, const DiscreteFunction &ubar ) {
            return self.assemble( ubar ).matrix().data();
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
        cls.def( "assemble", [] ( Scheme &self, const VirtualizedGridFunction< GridPart, RangeType > &ubar ) {
            return self.assemble( ubar ).matrix().data();
          }, pybind11::return_value_policy::reference_internal, "ubar"_a );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeAssemble ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeAssemble( cls, PriorityTag< 42 >() );
      }



      // registerSchemeGeneralCall
      // -------------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme & >()(
                     std::declval< const VirtualizedGridFunction< typename Scheme::GridPartType, typename Scheme::DiscreteFunctionSpaceType::RangeType > & >(),
                     std::declval< typename Scheme::DiscreteFunctionType & >()
                   ) ) >
      {
        typedef typename Scheme::DiscreteFunctionSpaceType::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        cls.def( "__call__", [] ( Scheme &self, const VirtualizedGridFunction< GridPart, RangeType > &arg, DiscreteFunction &dest ) { self( arg, dest ); } );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeGeneralCall ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeGeneralCall( cls, PriorityTag< 42 >() );
      }



      // registerSchemeModel
      // -------------------

      template< class Scheme, class... options >
      inline static auto registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 1 > )
        -> void_t< decltype( std::declval< Scheme >().model() ) >
      {
        cls.def_property_readonly( "model", &Scheme::model );
      }

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls, PriorityTag< 0 > )
      {}

      template< class Scheme, class... options >
      inline static void registerSchemeModel ( pybind11::class_< Scheme, options... > cls )
      {
        registerSchemeModel( cls, PriorityTag< 42 >() );
      }



      // registerScheme
      // --------------

      template< class Scheme, class... options >
      inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        using pybind11::operator""_a;

        registerSchemeConstructor( cls );

        cls.def( "_solve", [] ( Scheme &self, DiscreteFunction &solution ) {
            auto info = self.solve( solution );
            // needs pybind 1.9: return pybind11::dict("converged"_a=info.converged, "iterations"_a=info.nonlinearIterations, "linear_iterations"_a=info.linearIterations);
            return std::map<std::string,std::string> {
                {"converged",std::to_string(info.converged)},
                {"iterations",std::to_string(info.nonlinearIterations)},
                {"linear_iterations",std::to_string(info.linearIterations)}
              };
          } );
        cls.def( "__call__", [] ( Scheme &self, const DiscreteFunction &arg, DiscreteFunction &dest) { self( arg, dest ); } );
        registerSchemeGeneralCall( cls );

        cls.def_property_readonly( "dimRange", [] ( Scheme & ) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );
        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return detail::getSpace( self.cast< const Scheme & >(), self ); } );
        registerSchemeModel( cls );

        registerSchemeAssemble( cls );

        cls.def( "constraint", [] ( Scheme &self, DiscreteFunction &u) { self.constraint( u ); } );

        cls.def( "mark", [] ( Scheme &self, const DiscreteFunction &solution, double tolerance ) {
            double est = self.estimate( solution );
            return std::make_tuple( est, self.mark( tolerance ) );
          } );
      }

    } // namespace detail

    template< class Scheme, class... options >
    inline static void registerScheme ( pybind11::module module, pybind11::class_< Scheme, options... > cls )
    {
      detail::registerScheme( module, cls );
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SCHEME_HH
