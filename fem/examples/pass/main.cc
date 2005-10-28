// Global defines
#define DIM 3
#define DIM_OF_WORLD 3

#include "../../config.h"

// Include standard library headers
#include <iostream>
#include <string>

// Include fvcommon headers
#include "../../pass/dgpass.hh"

// Include Dune headers
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/io/file/grapedataio.hh>
#include "../../misc/inverseoperatorfactory.hh"
#include "../../misc/timestepping.hh"

using namespace Dune;
using namespace Adi;

template <class GridImp, 
          int dimRange, 
          int polOrd=0 >
struct Traits {
  typedef GridImp GridType;

  enum { dim = dimworld = GridType::dimension };

  typedef FunctionSpace < double, double , dim , dimRange > FuncSpace;
  typedef DofManager<GridType> DofManagerType;
  typedef LeafGridPart<GridType> GridPartType;
  typedef DiscontinuousFunctionSpace<FuncSpace,
                                     GridPartType,
                                     polOrd> FuncSpaceType ;
  typedef DFAdapt<FuncSpaceType> DiscreteFunctionType;
  typedef FieldVector<double, dimworld> NormalType;
  typedef FieldVector<double, dimRange> ValueType;
  typedef Identity<DiscreteFunctionType> IdentityType;

};

template <class TraitsImp>
class Fixture {
public:
  typedef TraitsImp Traits;
  typedef typename Traits::DomainType DomainType;

public:
  Fixture(std::string gridName, double epsilon, const DomainType& velo) :
    grid_(gridName),
    dm_(typename Traits::DofManagerFactory::getDofManager(grid_)),
    gridPart1_(grid_),
    gridPart2_(grid_),
    space1_(gridPart1_),
    space2_(gridPart2_),
    bc1_(), // * needs refinement
    bc2_(),
    problem1_(epsilon),
    problem2_(velo),
    pass0_(),
    pass1_(problem1_, pass0_, space1_, bc1_),
    pass2_(problem2_, pass1_, space2_, bc2_)
  {}

  typename Traits::Pass2Type& dgOperator() {
    return pass2_;
  }

private:
  typename Traits::GridType grid_;
  typename Traits::DofManagerType& dm_;
  typename Traits::GridPart1Type gridPart1_;
  typename Traits::GridPart2Type gridPart2_;
  typename Traits::Space1Type space1_;
  typename Traits::Space2Type space2_;
  typename Traits::BoundaryManager1Type bc1_;
  typename Traits::BoundaryManager2Type bc2_;
  typename Traits::Problem1Type problem1_;
  typename Traits::Problem2Type problem2_;
  typename Traits::Pass0Type pass0_;
  typename Traits::Pass1Type pass1_;
  typename Traits::Pass2Type pass2_;
};

template <class Loop>
void printData(double time, int id, Loop& loop) {

  std::cout << "Would print data" << std::endl;
}

int main() {
  typedef Traits<...> MyTraits;
  typedef Fixture<MyTraits> Fix;
  typedef MyTraits::Pass2Type OperatorType;
  typedef MyTraits::DomainType DomainType;
  typedef MyTraits::DiscreteFunctionType DiscreteFunction;

  // parameter section
  std::string gridFile("macro.small");
  DomainType velo(0.0); velo[0] = 1.0;
  double epsilon = 0.01;

  // problem creation
  Fix fix(gridFile, epsilon, velo);
  OperatorType& op = fix.dgOperator();

  // timestepping
  CGInverseOperatorFactory<DiscreteFunction> factory(1E-6, 1E-10, 100000, 0);
  ExplicitEuler<ProblemType> loop(problem, factory, dt);

  printData(loop.time(), 0, loop);
  loop.solveUntil(endTime);
  printData(loop.time(), 1, loop);
}
