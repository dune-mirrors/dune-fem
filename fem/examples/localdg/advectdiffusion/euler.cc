// Dune includes
#include <config.h>

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

#include <dune/common/misc.hh>
#include <dune/fem/common/boundary.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/quadraturerules.hh>

#include <iostream>
#include <string>

#include "advectdiff.hh"
#include "odesolver.hh"

using namespace Dune;
using namespace std;
#include "euler_mhd/eulermodel.hh"
#include "stuff.cc"
// Initial Data
class U0Smooth {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    if (arg[0]*arg[0] < 0.25) {
      // res[0] = cos(arg[0]*M_PI*2.)+2;
      res[0] = -8.*(arg[0]*arg[0]-0.25)+1.;
    }
    else {
      res[0] = 1.0;
    }
    res[1] = 1.5;
    res[2] = 0.;
    res[3] = 10.;
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
};
class U0VW {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    if (arg[0]<0.25) {
      res[0]=1.;
      res[1]=-1.;
      res[1]=0.;
      res[2]=0.;
      res[3]=1./(1.4-1.0);
    } else {
      res[0]=1.; 
      res[1]=1.;
      res[2]=0.;
      res[3]=1.0/(1.4-1.0);
    }
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
};
class U0Sod {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    if (arg[0]<0.25) {
      res[0]=1.;
      res[1]=0.;
      res[1]=0.;
      res[2]=0.;
      res[3]=1./(1.4-1.0);
    } else {
      res[0]=0.125;
      res[1]=0.;
      res[2]=0.;
      res[3]=0.1/(1.4-1.0);
    }
    res[1] *= res[0];
    res[2] *= res[0];
    res[3] += 0.5*res[1]*res[1]/res[0];
  }
};
int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  // Polynomial and ODE order
  enum {order=1,rksteps=2}; 
  const bool with_difftstep = false;
  // Grid:
  int N=40;                
  if (argc>1) 
    N = atoi(argv[1]);
  typedef SGrid<2, 2> GridType;
  SStruct s(N);
  GridType grid(s.n_, s.l_, s.h_);
  // CFL:
  double cfl;
  switch (order) {
  case 0: cfl=0.9;  break;
  case 1: cfl=0.2; break;
  case 2: cfl=0.1;  break;
  case 3: cfl=0.05;  break;
  case 4: cfl=0.09; break;
  }
  if (argc>2)
    cfl=atof(argv[2]);
  // *** Models
  typedef EulerModel<GridType> ModelType;
  ModelType euler(1.4,with_difftstep);
  // *** Fluxes 
  //typedef LLFFlux<ModelType> FluxType;
  typedef DWNumFlux<ModelType> FluxType;
  FluxType eulerflux(euler);
  // *** Operator typedefs
  // Space:
  //typedef DGAdvectionOperator<ModelType,LLFFlux,order> DgType;
  typedef DGLimitedAdvectionOperator<ModelType,DWNumFlux,order> DgType;
  typedef DuneODE::ExplTimeStepper<DgType> ODEType;
  DgType dg(grid,eulerflux);
  ODEType ode(dg,rksteps,cfl);
  // *** Initial data
  DgType::DestinationType U("U", dg.space());
  initialize<U0Smooth>(U);
  printSGrid(0,0,dg.space(),U);
  {
    DgType::DestinationType Utmp("Utmp", dg.space());
    dg(U,Utmp);
    printSGrid(0, 1, dg.space(), U);
  }
  
  // *** Time loop
  double t=0;
  while (t<1.) {
    t=ode.solve(U);
    cout << t << endl;
    ode.printGrid(1, U);
  }
  ode.printGrid(1, U);
} 




