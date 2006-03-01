// Dune includes
#include <config.h>

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/common/misc.hh>
#include <dune/fem/common/boundary.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/quadraturerules.hh>
#include <dune/io/file/grapedataio.hh>
#include <iostream>
#include <string>

using namespace Dune;
using namespace std;

#include "models.hh"
#include "stuff.cc"

int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  MPISTART
  // Polynomial and ODE order
  // Grid:
  GridType* grid=MacroGridParser().generate<GridType>(argv[1],MPI_COMM_WORLD);
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
  ModelType model(*grid,problem);
  // *** Fluxes 
  FluxType eulerflux(model);
  // *** Operator typedefs
  DgType dg(*grid,eulerflux);
  ODEType ode(dg,rksteps,cfl);
  // *** Initial data
  DgType::DestinationType U("U", dg.space());
  initialize(problem,U);
  // printSGrid(0,0,dg.space(),U);
  // *** Time loop
  double t=0;
  int n=0;
  double nextsave=0.;
  double savestep=0.05;
  double maxtime = 1.;
  while (t<maxtime) {
    t=ode.solve(U);
  }
  MPIEND
} 




