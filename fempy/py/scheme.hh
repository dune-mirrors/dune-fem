#ifndef DUNE_FEMPY_PY_SCHEME_HH
#define DUNE_FEMPY_PY_SCHEME_HH

#include <dune/fem/misc/l2norm.hh>

#include <dune/fempy/function/virtualizedgridfunction.hh>
#include <dune/fempy/py/common/numpyvector.hh>
#include <dune/fempy/py/discretefunction.hh>
#include <dune/fempy/pybind11/pybind11.h>

namespace Dune
{

  namespace FemPy
  {
    // registerScheme
    // -------------

    template< class Scheme >
    void registerScheme ( pybind11::module module )
    {
        typedef typename Scheme::DiscreteFunctionSpaceType Space;
        typedef typename Scheme::GridPartType GridPart;
        typedef typename Scheme::ModelType ModelType;
        typedef typename Scheme::DiscreteFunctionType DiscreteFunction;

        pybind11::class_< Scheme > cls( module, "Scheme" );

        cls.def( "__init__", [] ( Scheme &instance, Space &space, const ModelType& model, const std::string &prefix ) {
          new( &instance ) Scheme( space, model, prefix );
        }, pybind11::keep_alive< 1, 3 >(), pybind11::keep_alive< 1, 2 >() );
        cls.def("_solve", [] (Scheme &scheme, DiscreteFunction &solution) { return scheme.solve(solution, true); });
        cls.def("error", [] (Scheme &scheme, DiscreteFunction &solution)
        {
          const auto& gridExactSolution = scheme.exactSolution();
          Dune::Fem::L2Norm< GridPart > norm( solution.space().gridPart() );
          return norm.distance( solution, scheme.exactSolution() );
        } );
        cls.def_property_readonly( "name", &Scheme::name );
        cls.def_property_readonly( "dimRange", [](Scheme&) -> int { return DiscreteFunction::FunctionSpaceType::dimRange; } );
        cls.def_property_readonly( "space", &Scheme::space );
        /*
        cls.def("solution", [] (Scheme &scheme) -> DiscreteFunction& { return scheme.solution(); },
            pybind11::return_value_policy::reference_internal );
        */
    }
  }
}

#endif
