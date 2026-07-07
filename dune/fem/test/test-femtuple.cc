#include <config.h>

#include <dune/fem/function/tuplediscretefunction.hh>
#include <dune/fem/gridpart/adaptiveleafgridpart.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/space/discontinuousgalerkin/lagrange.hh>
#include <dune/grid/yaspgrid.hh>

using GridType = Dune::YaspGrid<2>;
using GridPartType = Dune::Fem::IntersectionAdaptiveLeafGridPart<GridType>;

using FunctionSpaceType = Dune::Fem::FunctionSpace<double, double, 2, 1>;
using SpaceType = Dune::Fem::LagrangeDiscontinuousGalerkinSpace<FunctionSpaceType, GridPartType, 1>;
using DiscreteFunctionType = Dune::Fem::AdaptiveDiscreteFunction<SpaceType>;
using TupleDiscreteFunctionType =
  Dune::Fem::TupleDiscreteFunction<DiscreteFunctionType, DiscreteFunctionType>;
using TupleSpaceType = typename TupleDiscreteFunctionType::DiscreteFunctionSpaceType;
using InverseOperatorType = Dune::Fem::KrylovInverseOperator<TupleDiscreteFunctionType>;

struct IdentityOperator
  : public Dune::Fem::Operator<TupleDiscreteFunctionType, TupleDiscreteFunctionType> {
  void operator()(const TupleDiscreteFunctionType& u, TupleDiscreteFunctionType& w) const {
    w.assign(u);
  }
};

int main(int argc, char** argv) try {
  Dune::Fem::MPIManager::initialize(argc, argv);

  GridType grid {{1.0, 1.0}, {4, 4}};
  GridPartType gridPart(grid);
  TupleSpaceType tupleSpace(gridPart);

  IdentityOperator linOp;
  TupleDiscreteFunctionType phi("solution", tupleSpace);
  TupleDiscreteFunctionType rhs("rhs", tupleSpace);
  phi.clear();
  rhs.clear();
  InverseOperatorType invOp(Dune::Fem::SolverParameter{});
  invOp.bind(linOp);
  invOp(rhs, phi);

  return 0;
} catch (const Dune::Exception& exception) {
  std::cerr << "Error: " << exception << std::endl;
  return 1;
}
