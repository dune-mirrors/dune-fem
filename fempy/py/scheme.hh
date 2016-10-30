#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fem/misc/l2norm.hh>

#include <dune/corepy/pybind11/pybind11.h>
#if HAVE_EIGEN
#include <dune/corepy/pybind11/eigen.h>
#endif

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
      // registerSchemeConstructor
      // -------------------------
      template< class Scheme, class Holder, class Alias >
      void registerSchemeConstructor ( pybind11::class_< Scheme, Holder, Alias > &cls, std::false_type )
      {}

      template< class Scheme, class Holder, class Alias >
      void registerSchemeConstructor ( pybind11::class_< Scheme, Holder, Alias > &cls, std::true_type )
      {
        using pybind11::operator""_a;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::ModelType ModelType;
        cls.def( "__init__", [] ( Scheme &self, Space &space, const ModelType &model ) {
            new (&self) Scheme( space, model );
          }, "space"_a, "model"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >() );
        cls.def( "__init__", [] ( Scheme &self, Space &space, const ModelType &model, const pybind11::dict &parameters ) {
            new (&self) Scheme( space, model, pyParameter( parameters, std::make_shared< std::string >() ) );
          }, "space"_a, "model"_a, "parameters"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >() );
      }

      template< class Scheme, class Holder, class Alias >
      void registerSchemeConstructor ( pybind11::class_< Scheme, Holder, Alias > &cls )
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::ModelType ModelType;
        registerSchemeConstructor( cls, std::is_constructible< Scheme, Space&, ModelType& >() );
      }

      // register assemble method if data method is available (and return value is registered)
      template <class Scheme,class Cls>
      auto registerSchemeAssemble(Cls &cls,int)
      // -> std::enable_if_t<Scheme::solver==SolverType::eigen,void>
      -> decltype(std::declval<typename Scheme::LinearOperatorType>().systemMatrix().matrix().data(),
          void())
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Space::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        cls.def("assemble", [](Scheme &scheme, const DiscreteFunction &ubar)
            { return scheme.assemble(ubar).systemMatrix().matrix().data(); },
            pybind11::return_value_policy::reference_internal
            );
        cls.def("assemble", [](Scheme &scheme, const VirtualizedGridFunction<GridPart,RangeType> &ubar)
            { return scheme.assemble(ubar).systemMatrix().matrix().data(); },
            pybind11::return_value_policy::reference_internal
            );
      }
      template <class Scheme,class Cls>
      void registerSchemeAssemble(Cls &cls,long)
      {}

      template <class Scheme, class Cls>
      auto registerSchemeGeneralCall( Cls &cls, int )
      -> decltype(std::declval<typename Scheme::DifferentiableOperatorType>().apply(
            std::declval<const VirtualizedGridFunction<typename Scheme::GridPartType,typename Scheme::DiscreteFunctionSpaceType::RangeType>&>(),
            std::declval<typename Scheme::DiscreteFunctionType&>() ),
          void())
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Space::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        cls.def("__call__", [] (Scheme &scheme,
                const VirtualizedGridFunction<GridPart,RangeType> &arg,
                DiscreteFunction &dest) { scheme(arg,dest); });
      }
      template< class Scheme, class Cls >
      auto registerSchemeGeneralCall( Cls &cls, long )
      {}

      template< class Scheme, class Cls >
      void registerScheme ( pybind11::module module, Cls &cls)
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        using pybind11::operator""_a;

        registerSchemeConstructor( cls );

        cls.def("_solve", [] (Scheme &scheme, DiscreteFunction &solution)
            { auto info = scheme.solve(solution);
              // needs pybind 1.9: return pybind11::dict("converged"_a=info.converged, "iterations"_a=info.nonlinearIterations, "linear_iterations"_a=info.linearIterations);
              return std::map<std::string,std::string> {
                  {"converged",std::to_string(info.converged)},
                  {"iterations",std::to_string(info.nonlinearIterations)},
                  {"linear_iterations",std::to_string(info.linearIterations)}
                };
            });
        cls.def("__call__", [] (Scheme &scheme, const DiscreteFunction &arg, DiscreteFunction &dest) { scheme(arg,dest); });
        registerSchemeGeneralCall<Scheme>(cls,0);

        cls.def_property_readonly( "dimRange", [](Scheme&) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );
        cls.def_property_readonly( "space", [] ( pybind11::object self ) { return detail::getSpace( self.cast< const Scheme & >(), self ); } );

        registerSchemeAssemble<Scheme>(cls,0);

        cls.def("constraint", [] (Scheme &scheme, DiscreteFunction &u) { scheme.constraint(u); });

        cls.def("mark", [] (Scheme &scheme, const DiscreteFunction &solution, double tolerance )
        {
          double est = scheme.estimate(solution);
          return std::make_tuple(est,scheme.mark(tolerance));
        });
      }
    }

    template< class Scheme, class Holder, class AliasType >
    void registerScheme ( pybind11::module module, pybind11::class_<Scheme, Holder, AliasType> &cls )
    {
      detail::registerScheme<Scheme>(module, cls);
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SCHEME_HH
