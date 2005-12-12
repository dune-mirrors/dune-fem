// Dune includes
#include "../../../config.h"

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
#include "scalarmodels.hh"
#include "stuff.cc"
// Initial Data
class U0 {
public:
  template <class DomainType, class RangeType>
  void evaluate(const DomainType& arg, RangeType& res) {
    double diffusion_ = 0.01;
    double t0 = 1.0;
    /*
    if (arg*arg < 0.25) {
      res = 1.0;
      res = cos(arg[0]*M_PI*2.)+1;
    }
    else {
      res = 0.0;
    }
    */

    res = 1./sqrt(4.*M_PI*diffusion_*t0)*exp(-(arg[0]+0.8-0.8*t0)*(arg[0]+0.8-0.8*t0)/(4.*diffusion_*t0));
  }
};
int main(int argc, char ** argv, char ** envp) {

  // *** Initialization


  // Polynomial and ODE order:
  enum {order=0,rksteps=1}; 
  const bool with_difftstep = true;
  

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


  // Diffusion parameter:
  double epsilon=0.01;
  if (argc>2)
    epsilon=atof(argv[2]);


  // *** Models
  // Advection-Diffusion Model
  typedef AdvectionDiffusionModel<GridType> AdvDiffType;

  // Define transport velocity
  AdvDiffType::DomainType velocity(0);
  velocity[0]=0.8;
  velocity[1]=0.;
  
  AdvDiffType advdiff(velocity,epsilon,with_difftstep);



  // *** Fluxes 
  // Advection Diffusion
  typedef LLFFlux<AdvDiffType> LLFAdvDiffType;
  LLFAdvDiffType llfadvdiff(advdiff);


  // *** DG and Time Integration Operators

  // Space (DG):
  typedef DGAdvectionDiffusionOperator<AdvDiffType,LLFFlux,order> DgAdvDiffType;
  
  // Time discretization:

  //      urspruenglich expl-rk methoden
  //typedef DuneODE::ExplRungeKutta<DgAdvDiffType> ODEAdvDiffType;

  //      dennis expl-rk methoden
  typedef DuneODE::ExplTimeStepper<DgAdvDiffType> ODEAdvDiffType;

  //      dennis impl-rk methoden
  //typedef DuneODE::ImplTimeStepper<DgAdvDiffType> ODEAdvDiffType;
  
  //      dennis semiimpl-rk methoden
  //typedef DuneODE::SemiImplTimeStepper<DgAdvType,DgDiffType> ODEAdvDiffType;
 


  // *** Construction...
  FieldVector<double,2> upwind;
  upwind[0]=1.;
  upwind[1]=1.;

  // Space:
  DgAdvDiffType dgadvdiff(grid,llfadvdiff,upwind);

  // Time:
  ODEAdvDiffType odeadvdiff(dgadvdiff,rksteps,cfl);


  // *** Initial data
  DgAdvDiffType::DestinationType U("U", dgadvdiff.space());
  initialize<U0>(U);



  // print initial data
  printSGrid(0, 1, dgadvdiff.space(), U);


  // *** Time loop
  double t=0,t1,dt,t_old;
  int n=0;t_old=0.;
  while (t<=1.) {
    n++;
    t=1.;
    if (odeadvdiff.time()<1.) {
      t1=odeadvdiff.solve(U);
      dgadvdiff.switchupwind();
      t=t1;
      dt= t - t_old; 
      t_old = t;
    }
    cout << n  << " : " <<  t << " : "  
	 << odeadvdiff.time() << " : " << dt << endl;
    odeadvdiff.printGrid(1, U);
  }

  // print final result
  odeadvdiff.printGrid(1, U);
} 
