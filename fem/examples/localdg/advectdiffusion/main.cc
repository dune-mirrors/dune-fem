// Dune includes
#include "../../../config.h"

#include "advectdiff.hh"

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

#include <dune/common/misc.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/common/boundary.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/quadraturerules.hh>

#include <iostream>
#include <string>

using namespace Dune;
using namespace std;

#include "scalarmodels.hh"

// ***********************
// ***********************************************
#include "stuff.cc"
class U0 {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    if (arg*arg < 0.25) {
      res = 1.0;
      res = cos(arg[0]*M_PI*2.)+1;
    }
    else {
      res = 0.0;
    }
    res*=100.; // (compensate some error in the initial projection)
  }
};
class TimeStepper : public TimeProvider {
 public:
  TimeStepper(double cfl) : 
    TimeProvider(0.), cfl_(cfl),
    savetime_(0.1), savestep_(1)
  {}
  template <class Operator>
  double solve(Operator &op,
	       typename Operator::DestinationType& U,
	       typename Operator::DestinationType& Upd) {
    resetTimeStepEstimate();
    op(U,Upd);
    double dt=cfl_*timeStepEstimate();
    //double dt=0.0001;
    Upd*=dt;
    U+=Upd;
    augmentTime(dt);
    return time();
  }
  template <class Operator>
  void printGrid(int nr, const typename Operator::SpaceType &space, 
		  const typename Operator::DestinationType& U) {
    if (time()>savetime_) {
      printSGrid(time(),savestep_*10+nr,space,U);
      ++savestep_;
      savetime_+=0.1;
    }
  }
 private:
  double cfl_;
  int savestep_;
  double savetime_;
};

int main() {
  enum {order=2};
  double cfl;
  switch (order) {
  case 0: cfl=1.0; break;
  case 1: cfl=0.1; break;
  case 2: cfl=0.05; break;
  case 3: cfl=0.1; break;
  case 4: cfl=0.09; break;
  }
  typedef SGrid<2, 2> GridType;
  SStruct s(200, 0.1);
  GridType grid(s.n_, s.l_, s.h_);

  // Model classes
  typedef AdvectionDiffusionModel<GridType> ModelType;
  typedef BurgersModel<GridType> BurgersType;
  ModelType::DomainType velocity(0.);
  velocity[0]=0.8;
  velocity[1]=0.;
  double epsilon=0.001;
  ModelType advdiff(velocity,0.);
  BurgersType burgers(epsilon);
  // Fluxes
  typedef UpwindFlux<ModelType> UpwindAdvDiffType;
  typedef LLFFlux<ModelType> LLFAdvDiffType;
  typedef LLFFlux<BurgersType> LLFBurgers;
  UpwindAdvDiffType upwindadvdiff(advdiff);
  LLFAdvDiffType llfadvdiff(advdiff);
  LLFBurgers llfburgers(burgers);
  // ODE Solvers
  typedef TimeStepper ODEType;
  ODEType ode(cfl*0.9);
  ODEType odeLLF(cfl*0.9);
  ODEType odeburgers(cfl*0.9);
  // Operators
  typedef DGAdvectionDiffusionOperator<BurgersType,LLFFlux,order> DgTypeBurgers;
  typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgTypeLLF;
  typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
  FieldVector<double,2> upwind;
  upwind[0]=1.;
  upwind[1]=1.;
  DgType dg(grid,ode,upwindadvdiff,upwind);
  DgTypeLLF dgLLF(grid,odeLLF,llfadvdiff,upwind);
  DgTypeBurgers dgBurgers(grid,odeburgers,llfburgers,upwind);

  // Grid and Data...
  DgType::DestinationType U("U", dg.space());
  DgType::DestinationType ULLF("ULLF", dg.space());
  DgType::DestinationType UBurgers("UBurgers", dg.space());
  DgType::DestinationType Upd("Upd", dg.space());
  initialize<U0>(U);
  initialize<U0>(ULLF);
  initialize<U0>(UBurgers);
  printSGrid(0, 1, dg.space(), U);
  printSGrid(0, 2, dg.space(), ULLF);
  printSGrid(0, 3, dg.space(), UBurgers);
  double t=0,t1,t2,t3;
  while (t<1.) {
    if (ode.time()<1.)
      t1=ode.solve(dg,U,Upd);
    if (odeLLF.time()<1.)
      t2=odeLLF.solve(dgLLF,ULLF,Upd);
    if (odeburgers.time()<1.)
      t3=odeburgers.solve(dgBurgers,UBurgers,Upd);
    t=(t1<t2)?t1:t2;
    t=(t<t3)?t:t3;
    cout << t << " : "  
	 << ode.time() << " " << odeLLF.time() << " " << odeburgers.time() << endl;
    ode.printGrid<DgType>(1,dg.space(), U);
    odeLLF.printGrid<DgType>(2,dg.space(), ULLF);
    odeburgers.printGrid<DgType>(3,dg.space(), UBurgers);
  }
  ode.printGrid<DgType>(1,dg.space(), U);
  odeLLF.printGrid<DgType>(2,dg.space(), ULLF);
  odeburgers.printGrid<DgType>(3,dg.space(), UBurgers);
}
