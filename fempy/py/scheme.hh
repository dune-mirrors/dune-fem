#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fem/misc/l2norm.hh>
#include <dune/fem/schemes/solver.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/parameter.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/corepy/pybind11/pybind11.h>
#if HAVE_EIGEN
  #include <dune/corepy/pybind11/eigen.h>
#endif

namespace Dune
{

  namespace FemPy
  {
    // registerScheme
    // -------------

    namespace detail
    {
      // register assemble method if data method is available (and return value is registered)
      template <class Scheme,class Cls>
      auto registerSchemeAssemble(Cls &cls,int)
      -> std::enable_if_t<Scheme::solver==SolverType::eigen,void>
      {
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Space::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;

        cls.def("assemble", [](Scheme &scheme, const DiscreteFunction &ubar)
            { return scheme.assemble(ubar).systemMatrix().matrix().data(); } );
            // pybind11::keep_alive<0,1>());
        cls.def("assemble", [](Scheme &scheme, const VirtualizedGridFunction<GridPart,RangeType> &ubar)
            { return scheme.assemble(ubar).systemMatrix().matrix().data(); } );
            // pybind11::keep_alive<0,1>());
      }
      template <class Scheme,class Cls>
      void registerSchemeAssemble(Cls &cls,long)
      {}

      template< class Scheme, class Cls >
      void registerScheme ( pybind11::module module, Cls &cls, std::true_type)
      {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Space::RangeType RangeType;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::ModelType ModelType;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        using pybind11::operator""_a;

        cls.def( "__init__", [] ( Scheme &self, Space &space, const ModelType &model, std::string name ) {
            new (&self) Scheme( space, model, std::move( name ) );
          }, "space"_a, "model"_a, "name"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >() );

        cls.def( "__init__", [] ( Scheme &self, Space &space, const ModelType &model, std::string name, const pybind11::dict &parameters ) {
            new (&self) Scheme( space, model, std::move( name ), pyParameter( parameters, std::make_shared< std::string >() ) );
          }, "space"_a, "model"_a, "name"_a, "parameters"_a, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >() );

        cls.def("_prepare", [] (Scheme &scheme, const DiscreteFunction &add) { scheme.prepare(add); });
        cls.def("_prepare", [] (Scheme &scheme) { scheme.prepare(); });
        cls.def("_solve", [] (Scheme &scheme, DiscreteFunction &solution,bool assemble) { scheme.solve(solution, assemble); });
        cls.def("__call__", [] (Scheme &scheme, const DiscreteFunction &arg, DiscreteFunction &dest) { scheme(arg,dest); });
        cls.def("__call__", [] (Scheme &scheme, const VirtualizedGridFunction<GridPart,RangeType> &arg, DiscreteFunction &dest) { scheme(arg,dest); });

        cls.def("error", [] (Scheme &scheme, DiscreteFunction &solution)
        {
          Dune::Fem::L2Norm< GridPart > norm( solution.space().gridPart() );
          return norm.distance( solution, scheme.exactSolution() );
        } );
        cls.def("mark", [] (Scheme &scheme, const DiscreteFunction &solution, double tolerance )
        {
          double est = scheme.estimate(solution);
          return std::make_tuple(est,scheme.mark(tolerance));
        });
        cls.def_property_readonly( "name", &Scheme::name );
        cls.def_property_readonly( "dimRange", [](Scheme&) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );

        cls.def_property_readonly( "space", [] ( const Scheme &self ) {
            const Space &space = self.space();
            pybind11::handle hSpace = pybind11::detail::get_object_handle( &space, pybind11::detail::get_type_info( typeid( Space ) ) );
            return pybind11::object( hSpace, true );
          } );

        registerSchemeAssemble<Scheme>(cls,0);

        cls.def("constraint", [] (Scheme &scheme, DiscreteFunction &u) { scheme.constraint(u); });
      }
    }

    template< class Scheme, class Holder, class AliasType >
    void registerScheme ( pybind11::module module, pybind11::class_<Scheme, Holder, AliasType> &cls )
    {
      typedef typename Scheme::DiscreteFunctionSpaceType Space;
      typedef typename Scheme::ModelType ModelType;
      detail::registerScheme<Scheme>(module, cls, std::is_constructible<Scheme, Space&, ModelType&, std::string&>());
    }

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_PY_SCHEME_HH
