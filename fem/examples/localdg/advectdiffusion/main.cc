// #define DUNE_DEPRECATED
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
#include <dune/fem/l2projection.hh>

#include <iostream>
#include <string>

using namespace Dune;
using namespace std;

#include "scalarmodels.hh"
template <class GridType>
double BurgersModel<GridType>::epsilon;
template <class GridType>
typename AdvectionDiffusionModel<GridType>::DomainType 
       AdvectionDiffusionModel<GridType>::velocity;
template <class GridType>
double AdvectionDiffusionModel<GridType>::epsilon;

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
  TimeStepper(double cfl) : TimeProvider(0.), cfl_(cfl) {}
  template <class Operator>
  double solve(Operator &op,
	       typename Operator::DestinationType& U,
	       typename Operator::DestinationType& Upd) {
    resetTimeStepEstimate();
    op(U,Upd);
    // double dt=cfl_*timeStepEstimate();
    double dt=0.0001;
    Upd*=dt;
    U+=Upd;
    augmentTime(dt);
    return time();
  }
 private:
  double cfl_;
};
template <class Operator>
double solve(Operator &op,
	     typename Operator::DestinationType& U,
	     typename Operator::DestinationType& Upd) {
  Upd.clear();
  op(U,Upd);
  double dt=0.001;
  Upd*=dt;
  U+=Upd;
  return dt;
}

int main() {
  enum {order=0};
  typedef SGrid<2, 2> GridType;
  typedef AdvectionDiffusionModel<GridType> ModelType;
  typedef BurgersModel<GridType> BurgersType;
  ModelType::velocity[0]=0.8;
  ModelType::velocity[1]=0.;
  ModelType::epsilon=0.01;
  BurgersType::epsilon=0.01;
  
  typedef DGAdvectionDiffusionOperator<BurgersType,LLFFlux,order> DgTypeBurgers;
  typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgTypeLLF;
  typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
  typedef TimeStepper ODEType;

  SStruct s(200, 0.1);
  GridType grid(s.n_, s.l_, s.h_);
  FieldVector<double,2> upwind;
  upwind[0]=1.;
  upwind[1]=1.;

  ODEType ode(0.1);
  ODEType odeLLF(0.1);
  ODEType odeburgers(0.1);

  DgType dg(grid,ode,upwind);
  DgTypeLLF dgLLF(grid,odeLLF,upwind);
  DgTypeBurgers dgBurgers(grid,odeburgers,upwind);

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
  double save=0.1;
  int step=1;
  double t=0;
  while (t<1.) {
    /*
    t+=solve(dg,U,Upd);
    solve(dgLLF,ULLF,Upd);
    solve(dgBurgers,UBurgers,Upd);
    */
    // t=ode.solve(dg,U,Upd);
    t+=solve(dg,U,Upd);
    cout << t << endl;
    if (t>save) {
      printSGrid(0, 10*step+1, dg.space(), U);
      printSGrid(0, 10*step+2, dg.space(), ULLF);
      printSGrid(0, 10*step+3, dg.space(), UBurgers);
      save+=0.1;
      ++step;
    }
  }
  printSGrid(0, 10*step+1, dg.space(), U);
  printSGrid(0, 10*step+2, dg.space(), ULLF);
  printSGrid(0, 10*step+3, dg.space(), UBurgers);
}
