// Dune includes
#include <config.h>

#include <dune/fem/pass/utility.hh>

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/common/misc.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <iostream>
#include <string>
using namespace std;

#include "models.hh"
#include "advectdiff.hh"

#include <dune/fem/solver/odesolver.hh>

using namespace Dune;
//#include "scalarmodels.hh"			
#include "scalarmodels.hh"
#include "stuff.cc"
// Initial Data
int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  // Polynomial and ODE order
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
  // Initial Data and Problem Description
  U0<GridType> problem1(epsilon,with_difftstep);
  U0<GridType> problem2(0,false);
  // *** Models
  // Advection-Diffusion Model
  typedef AdvectionDiffusionModel<GridType,U0<GridType> > AdvDiffType;
  // Advection-Diffusion
  AdvDiffType advdiff(grid,problem1);
  // Diffusion
  AdvDiffType diffeqn(grid,problem1);
  // Advection
  AdvDiffType adveqn(grid,problem2);  
  // Burgers Model
  typedef BurgersModel<GridType,U0<GridType> > BurgersType;
  BurgersType burgers(grid,problem1);
  typedef BurgersModel<GridType,U0<GridType> > BurgersAdvType;
  BurgersType burgersadv(grid,problem2);
  // *** Fluxes 
  // Advection Diffusion
  typedef UpwindFlux<AdvDiffType> UpwindAdvDiffType;
  typedef LLFFlux<AdvDiffType> LLFAdvDiffType;
  UpwindAdvDiffType llfadvdiff(advdiff);
  LLFAdvDiffType llfdiffeqn(diffeqn);
  UpwindAdvDiffType upwindadveqn(adveqn);
  // Burgers Model
  typedef LLFFlux<BurgersType> LLFBurgers;
  LLFBurgers llfburgers(burgers);
  typedef LLFFlux<BurgersAdvType> LLFBurgersAdv;
  LLFBurgersAdv llfburgersadv(burgersadv);
  // *** Operator typedefs
  // Space:
  typedef DGAdvectionDiffusionOperator<AdvDiffType,UpwindFlux,order> DgAdvDiffType;
  typedef DGDiffusionOperator<AdvDiffType,LLFFlux,order> DgDiffType;
  typedef DGAdvectionOperator<AdvDiffType,UpwindFlux,order> DgAdvType;
  typedef DGAdvectionDiffusionOperator<BurgersType,LLFFlux,order> DgBurgersType;
  typedef DGAdvectionOperator<BurgersAdvType,LLFFlux,order> DgBurgersAdvType;
  // Time:
  //      urspruenglich expl-rk methoden
  /*
  typedef DuneODE::ExplRungeKutta<DgAdvDiffType> ODEAdvDiffType;
  typedef DuneODE::ExplRungeKutta<DgDiffType> ODEDiffType;
  typedef DuneODE::ExplRungeKutta<DgAdvType> ODEAdvType;
  typedef DuneODE::ExplRungeKutta<DgBurgersType> ODEBurgersType;
  */
  //      dennis expl-rk methoden
  typedef DuneODE::ExplTimeStepper<DgAdvDiffType> ODEAdvDiffType;
  typedef DuneODE::ExplTimeStepper<DgDiffType> ODEDiffType;
  typedef DuneODE::ExplTimeStepper<DgAdvType> ODEAdvType;
  typedef DuneODE::ExplTimeStepper<DgBurgersType> ODEBurgersType;
  //      dennis impl-rk methoden
  /*
  typedef DuneODE::ImplTimeStepper<DgAdvDiffType> ODEAdvDiffType;
  typedef DuneODE::ImplTimeStepper<DgDiffType> ODEDiffType;
  typedef DuneODE::ImplTimeStepper<DgAdvType> ODEAdvType;
  typedef DuneODE::ImplTimeStepper<DgBurgersType> ODEBurgersType;
  */
  //      dennis semiimpl-rk methoden
  /*
  typedef DuneODE::SemiImplTimeStepper<DgAdvType,DgDiffType> ODEAdvDiffType;
  typedef DuneODE::ImplTimeStepper<DgDiffType> ODEDiffType;
  typedef DuneODE::ImplTimeStepper<DgAdvType> ODEAdvType;
  typedef DuneODE::SemiImplTimeStepper<DgBurgersAdvType,DgDiffType> ODEBurgersType;
  */
  // *** Construction...
  FieldVector<double,2> upwind;
  upwind[0]=1.;
  upwind[1]=1.;
  // Space:
  DgAdvDiffType dgadvdiff(grid,llfadvdiff,upwind);
  DgDiffType dgdiffeqn(grid,llfdiffeqn,upwind);
  DgAdvType dgadveqn(grid,upwindadveqn);
  DgAdvType dgadv1eqn(grid,upwindadveqn);            // for semi-impl
  DgBurgersType dgburgers(grid,llfburgers,upwind); 
  DgBurgersAdvType dgburgersadv(grid,llfburgersadv); 
  // Time:
  ODEAdvDiffType odeadvdiff(dgadvdiff,rksteps,cfl);
  // ODEAdvDiffType odeadvdiff(dgadv1eqn,dgdiffeqn,rksteps,cfl); // semi-impl
  ODEDiffType odediffeqn(dgdiffeqn,rksteps,cfl);
  ODEAdvType odeadveqn(dgadveqn,rksteps,cfl);
  ODEBurgersType odeburgers(dgburgers,rksteps,cfl);
  // ODEBurgersType odeburgers(dgburgersadv,dgdiffeqn,rksteps,cfl); // semi-impl
  // *** Initial data
  DgAdvDiffType::DestinationType U("U", dgadvdiff.space());
  DgDiffType::DestinationType UDiff("DDiff", dgdiffeqn.space());
  DgAdvType::DestinationType UAdv("UAdv", dgadveqn.space());
  DgBurgersType::DestinationType UBurgers("UBurgers", dgburgers.space());
  // Initialize data
  initialize(problem1,U);
  initialize(problem1,UDiff);
  initialize(problem1,UAdv);
  initialize(problem1,UBurgers);

  printSGrid(0, 1, dgadvdiff.space(), U);
  printSGrid(0, 2, dgdiffeqn.space(), UDiff);
  printSGrid(0, 3, dgadveqn.space(), UAdv);
  printSGrid(0, 4, dgburgers.space(), UBurgers);
  // *** Time loop
  double t=0,t1,t2,t3,t4;
  while (t<1.) {
    t=1.;
    if (odeadvdiff.time()<1.) {
      t1=odeadvdiff.solve(U);
      dgadvdiff.switchupwind();
      t=t1;
    }
    if (odediffeqn.time()<1.) {
      t2=odediffeqn.solve(UDiff);
      dgdiffeqn.switchupwind();
      t=(t<t2)?t:t2;
    } 
    if (odeadveqn.time()<1.) {
      t3=odeadveqn.solve(UAdv);
      t=(t<t3)?t:t3;
    }
    if (odeburgers.time()<1.) {
      t4=odeburgers.solve(UBurgers);
      dgburgers.switchupwind();
      t=(t<t4)?t:t4;
    }
    cout << t << " : "  
	 << odeadvdiff.time() << " " << odediffeqn.time() << " " 
	 << odeadveqn.time() << " " << odeburgers.time() << endl;
    odeadvdiff.printGrid(1, U);
    odediffeqn.printGrid(2, UDiff);
    odeadveqn.printGrid(3, UAdv);
    odeburgers.printGrid(4, UBurgers);
  }
  odeadvdiff.printGrid(1, U);
  odediffeqn.printGrid(2, UDiff);
  odeadveqn.printGrid(3, UAdv);
  odeburgers.printGrid(4, UBurgers);
} 
