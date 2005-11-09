// Global defines
#undef DIM
#undef DIM_OF_WORLD

#define DIM 2
#define DIM_OF_WORLD 2

#include "../../config.h"

// Include standard library headers
#include <iostream>
#include <string>

// Include fvcommon headers
#include "../../pass/dgpass.hh"
#include "problem.hh"

// Include Dune headers
#include <dune/common/misc.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/common/boundary.hh>
//#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/albertagrid.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/io/file/grapedataio.hh>
#include <dune/fem/l2projection.hh>
#include "../../misc/inverseoperatorfactory.hh"
#include "../../misc/timenew.hh"
#include "../../misc/identity.hh"

// Dx includes
#include "../../visual/dx/dxdata.hh"

using namespace Dune;
//using namespace Adi;

template <class GridImp, 
          int polOrd >
struct Traits {
  typedef GridImp GridType;

  enum { dim = GridType::dimension };
  enum { dimRange = 1 };

  typedef FunctionSpace < double, double , dim , dimRange> FuncSpace;
  typedef DofManager<GridType> DofManagerType;
  typedef DofManagerFactory<DofManagerType> DofManagerFactory;
  //typedef LeafGridPart<GridType> GridPartType;
  typedef DefaultGridIndexSet<GridType, LevelIndex> IndexSetType;
  typedef DefaultGridPart<GridType, IndexSetType> GridPartType;
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
  typedef IndexSetType GridIndexSet1;
  typedef DiscontinuousGalerkinSpace<
    SingleFuncSpace1, GridPartType, polOrd> SingleSpace1Type;
  typedef CombinedSpace<SingleSpace1Type, dim> Space1Type;
  typedef TransportDiffusionProblem1 Problem1Type;
  typedef LocalDGPass<Problem1Type, Pass0Type> Pass1Type;

  // * types for pass2
  typedef FunctionSpace <double, double , dim , dimRange> FuncSpace2;
  typedef GridPartType GridPart2Type;
  typedef IndexSetType GridIndexSet2;
  typedef DiscontinuousGalerkinSpace<
    FuncSpace2, GridPartType, polOrd> Space2Type;
  typedef TransportDiffusionProblem2 Problem2Type;
  typedef LocalDGPass<Problem2Type, Pass1Type> Pass2Type;
};

struct SStruct {
  SStruct(int n1, int n2, double lx, double ly, double hx, double hy) 
  {
    n_[0] = n1;
    n_[1] = n2;
    l_[0] = lx;
    l_[1] = ly;
    h_[0] = hx;
    h_[1] = hy;
  }

  SStruct(int n, double h) {
    n_[0] = n;
    n_[1] = 1;
    l_[0] = -1.0;
    l_[1] = -h/2.0;
    h_[0] = 1.0;
    h_[1] = h/2.0;
  }

  int n_[2];
  double l_[2];
  double h_[2];
};

template <class TraitsImp>
class Fixture {
public:
  typedef TraitsImp Traits;
  typedef typename Traits::DomainType DomainType;

public:
  Fixture(std::string gridName, double epsilon, const DomainType& velo, SStruct& s) :
    //    grid_(gridName),
    grid_(s.n_, s.l_, s.h_),
    dm_(Traits::DofManagerFactory::getDofManager(grid_)),
    iset1_(grid_, grid_.maxLevel()),
    iset2_(grid_, grid_.maxLevel()),
    gridPart1_(grid_, iset1_),
    gridPart2_(grid_, iset2_),
    singleSpace1_(gridPart1_),
    space1_(singleSpace1_),
    space2_(gridPart2_),
    problem1_(),
    problem2_(velo, epsilon),
    pass0_(),
    pass1_(problem1_, pass0_, space1_),
    pass2_(problem2_, pass1_, space2_)
  {}

  typename Traits::GridType& grid() {
    return grid_;
  }

  typename Traits::DofManagerType& dm() {
    return dm_;
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
  typename Traits::GridIndexSet1 iset1_;
  typename Traits::GridIndexSet2 iset2_;
  typename Traits::GridPart1Type gridPart1_;
  typename Traits::GridPart2Type gridPart2_;
  typename Traits::SingleSpace1Type singleSpace1_;
  typename Traits::Space1Type space1_;
  typename Traits::Space2Type space2_;
  typename Traits::Problem1Type problem1_;
  typename Traits::Problem2Type problem2_;
  typename Traits::Pass0Type pass0_;
  typename Traits::Pass1Type pass1_;
  typename Traits::Pass2Type pass2_;
};

/*
template <class Loop, class GridType>
void printData(double time, int timestep, GridType& grid, Loop& loop) 
{
  GrapeDataIO<GridType> dataio;
  dataio.writeGrid(grid, xdr, "grid", time, timestep);
  dataio.writeData(loop.solution(), xdr, "vec", timestep);
}
*/

template <class Loop, class SpaceType>
void printDX(double time, int timestep, SpaceType& space, Loop& loop) 
{  
  std::ostringstream stream, filestream;
  filestream << "sol" << timestep;
  DXWriter<SpaceType, false> dx(space, filestream.str());
  dx.write(loop.solution(), "data");
}

template <class Geometry>
void midPoint(const Geometry& geo, FieldVector<double, 2>& result) 
{
  result *= 0.0;
  for (int i = 0; i < geo.corners(); ++i) {
    result += geo[i];
  }

  result /= static_cast<double>(geo.corners());
}

template <class Loop, class SpaceType>
void printSGrid(double time, int timestep, SpaceType& space, Loop& loop)
{
  typedef typename SpaceType::IteratorType Iterator;
  typedef typename Loop::DiscreteFunctionType DiscreteFunctionType;
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;
    
  std::ostringstream filestream;
  filestream << "sgrid" << timestep;
  
  std::ofstream ofs(filestream.str().c_str(), std::ios::out);

  FieldVector<double, 2> mid(0.0);
  FieldVector<double, 2> localMid(0.5);

  FieldVector<double, 1> result;

  DiscreteFunctionType& sol = loop.solution();

  Iterator endit = space.end();
  for (Iterator it = space.begin(); it != endit; ++it) {
    midPoint(it->geometry(), mid);
    if (mid[1] < 0.1 && mid[1] > -0.1) {
      LocalFunctionType lf = sol.localFunction(*it);
      lf.evaluateLocal(*it, localMid, result);
      ofs << mid[0] << " " << result[0] << "\n";
    }
  }
  ofs << std::endl;
}

class StupidFunction {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    if (arg*arg < 0.25) {
      res = 1.0;
    } 
    else {
      res = 0.0;
    }
  }
};


template <class DFType>
void initialize(DFType& df) 
{
  //- Typedefs and enums
  typedef typename DFType::Traits::DiscreteFunctionSpaceType SpaceType;
  typedef typename SpaceType::Traits::IteratorType Iterator;
  typedef typename DFType::Traits::LocalFunctionType LocalFunction;
  typedef typename SpaceType::Traits::GridType Grid;
  typedef typename Grid::template Codim<0>::Entity::Geometry Geometry;

  typedef typename SpaceType::Traits::RangeType RangeType;
  typedef typename SpaceType::Traits::DomainType DomainType;

  enum { dim = Grid::dimension };

  typedef FieldVector<double, dim> Coordinate;

  //- Local classes

  //- Actual method
  StupidFunction f;  
  L2Projection<DFType, StupidFunction, 2>::project(f, df);

  typedef typename DFType::DofIteratorType DofIterator;
  for (DofIterator it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
  }
}

template <class DFType>
void printIt(DFType& df) 
{
  typedef typename DFType::DofIteratorType DofIt;

  std::cout << "print it\n";
  for (DofIt it = df.dbegin(); it != df.dend(); ++it) {
    std::cout << *it << std::endl;
  }

}

int main() 
{

  typedef SGrid<2, 2> MyGrid;

  //  typedef AlbertaGrid<2, 2> MyGrid;

  typedef Traits<MyGrid, 3> MyTraits;
  typedef Fixture<MyTraits> Fix;
  typedef MyTraits::Pass2Type SpaceOperatorType;
  typedef MyTraits::DomainType DomainType;
  typedef MyTraits::DiscreteFunctionType DiscreteFunction;
  typedef MyTraits::DestinationType DestinationType;
  typedef MyTraits::MappingType MappingType;
  typedef MyTraits::IdentityType TimeOperatorType;
  typedef MyTraits::GridType GridType;

  // parameter section
  /*
  std::string gridFile("macro.small");
  DomainType velo(0.0); velo[0] = 1.0;
  double epsilon = 0.01;
  double dt = 0.1;
  double endTime = 1.0;
  int globalRefineCount = 2;

  */
  std::string gridFile;
  DomainType velo(0.0);
  double epsilon;
  double dt;
  double endTime;
  int globalRefineCount;
  
  double veloX;
  std::string paramfile("parameter");
  readParameter(paramfile,"gridFile", gridFile);
  readParameter(paramfile,"epsilon", epsilon);
  readParameter(paramfile,"velo", veloX);
  velo[0] = veloX;
  readParameter(paramfile,"dt", dt);
  readParameter(paramfile,"endTime", endTime);
  readParameter(paramfile,"globalRefineCount", globalRefineCount);

  // problem creation
  SStruct s(100, 0.1);
  Fix fix(gridFile, epsilon, velo, s);
  GridType& grid = fix.grid();

  grid.globalRefine(globalRefineCount);
  fix.dm().resize();
  fix.dm().dofCompress();

  //grid.globalRefine(globalRefineCount);
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
  op.setTimeProvider(loop);

  //printData(loop.time(), 0, fix.grid(), loop);
  printSGrid(loop.time(), 0, fix.space(), loop);
  printDX(loop.time(), 0, fix.space(), loop);
  
  loop.solve();
  //printData(loop.time(), 1, fix.grid(), loop);
  printSGrid(loop.time(), 1, fix.space(), loop);
  printDX(loop.time(), 1, fix.space(), loop);

  loop.solveUntil(endTime);
  //printData(loop.time(), 2, fix.grid(), loop);
  printSGrid(loop.time(), 2, fix.space(), loop);
  printDX(loop.time(), 2, fix.space(), loop);

}
