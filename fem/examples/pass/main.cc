// Global defines
#define DIM 3
#define DIM_OF_WORLD 3

#include "../../config.h"

// Include standard library headers
#include <iostream>
#include <string>

// Include fvcommon headers
#include "../../pass/dgpass.hh"
#include "problem.hh"

// Include Dune headers
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/common/boundary.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/io/file/grapedataio.hh>
#include "../../misc/inverseoperatorfactory.hh"
#include "../../misc/timenew.hh"
#include "../../misc/identity.hh"

using namespace Dune;
//using namespace Adi;

template <class GridImp, 
          int polOrd=0 >
struct Traits {
  typedef GridImp GridType;

  enum { dim = GridType::dimension };
  enum { dimRange = 1 };

  typedef FunctionSpace < double, double , dim , dimRange> FuncSpace;
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactory;
  typedef LeafGridPart<GridType> GridPartType;
  typedef DiscontinuousGalerkinSpace<FuncSpace,
                                     GridPartType,
                                     polOrd> FuncSpaceType ;
  typedef AdaptiveDiscreteFunction<FuncSpaceType> DiscreteFunctionType;
  typedef FieldVector<double, dim> DomainType;
  typedef FieldVector<double, dimRange> RangeType;
  typedef DiscreteFunctionType DestinationType;
  typedef Identity<DiscreteFunctionType> IdentityType;
  typedef StartPass<DiscreteFunctionType> Pass0Type;
  typedef Mapping<
    double, double, DestinationType, DestinationType> MappingType;

  // * types for pass1
  typedef FunctionSpace<double, double , dim , dim> FuncSpace1;
  typedef FunctionSpace<double, double, dim, 1> SingleFuncSpace1;
  typedef GridPartType GridPart1Type;
  typedef DiscontinuousGalerkinSpace<
    SingleFuncSpace1, GridPartType, polOrd> SingleSpace1Type;
  typedef CombinedSpace<SingleSpace1Type, dim> Space1Type;
  typedef BoundaryManager<Space1Type> BoundaryManager1Type;
  typedef TransportDiffusionProblem1 Problem1Type;
  typedef LocalDGPass<Problem1Type, Pass0Type> Pass1Type;

  // * types for pass2
  typedef FunctionSpace <double, double , dim , dimRange> FuncSpace2;
  typedef GridPartType GridPart2Type;
  typedef DiscontinuousGalerkinSpace<
    FuncSpace2, GridPartType, polOrd> Space2Type;
  typedef BoundaryManager<Space2Type> BoundaryManager2Type;
  typedef TransportDiffusionProblem2 Problem2Type;
  typedef LocalDGPass<Problem2Type, Pass1Type> Pass2Type;
};

template <class TraitsImp>
class Fixture {
public:
  typedef TraitsImp Traits;
  typedef typename Traits::DomainType DomainType;

public:
  Fixture(std::string gridName, double epsilon, const DomainType& velo) :
    grid_(gridName),
    dm_(Traits::DofManagerFactory::getDofManager(grid_)),
    gridPart1_(grid_),
    gridPart2_(grid_),
    singleSpace1_(gridPart1_),
    space1_(singleSpace1_),
    space2_(gridPart2_),
    bc1_(), // * needs refinement
    bc2_(),
    problem1_(epsilon),
    problem2_(velo),
    pass0_(),
    pass1_(problem1_, pass0_, space1_, bc1_),
    pass2_(problem2_, pass1_, space2_, bc2_)
  {}

  typename Traits::GridType& grid() {
    return grid_;
  }

  typename Traits::Space2Type& space() {
    return space2_;
  }

  typename Traits::Pass2Type& dgOperator() {
    return pass2_;
  }

private:
  typename Traits::GridType grid_;
  typename Traits::DofManagerType& dm_;
  typename Traits::GridPart1Type gridPart1_;
  typename Traits::GridPart2Type gridPart2_;
  typename Traits::SingleSpace1Type singleSpace1_;
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

template <class Loop, class GridType>
void printData(double time, int timestep, GridType& grid, Loop& loop) 
{
  GrapeDataIO<GridType> dataio;
  dataio.writeGrid(grid, xdr, "grid", time, timestep);
  dataio.writeData(loop.solution(), xdr, "vec", timestep);
}

template <class DFType>
void initialize(DFType& df) 
{
  typedef typename DFType::DofIteratorType Iterator;
  for (Iterator it = df.dbegin(); it != df.dend(); ++it) {
    *it = 1.0;
  }
}

int main() 
{
  typedef Traits<ALU3dGrid<3, 3, tetra>, 0> MyTraits;
  typedef Fixture<MyTraits> Fix;
  typedef MyTraits::Pass2Type SpaceOperatorType;
  typedef MyTraits::DomainType DomainType;
  typedef MyTraits::DiscreteFunctionType DiscreteFunction;
  typedef MyTraits::DestinationType DestinationType;
  typedef MyTraits::MappingType MappingType;
  typedef MyTraits::IdentityType TimeOperatorType;

  // parameter section
  std::string gridFile("macro.small");
  DomainType velo(0.0); velo[0] = 1.0;
  double epsilon = 0.01;
  double dt = 0.1;
  double endTime = 1.0;

  // problem creation
  Fix fix(gridFile, epsilon, velo);
  SpaceOperatorType& op = fix.dgOperator();
  TimeOperatorType top;

  // intial data
  DestinationType initial("initial", fix.space());
  initialize(initial);

  // timestepping
  IdentitySolverFactory<DiscreteFunction> factory;
  //CGInverseOperatorFactory<DiscreteFunction> factory(1E-6, 1E-10, 100000, 0);
  ExplicitEuler<
    DiscreteFunction, MappingType> loop(initial, top, op, factory, dt);

  printData(loop.time(), 0, fix.grid(), loop);
  loop.solveUntil(endTime);
  printData(loop.time(), 1, fix.grid(), loop);
}
