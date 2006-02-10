// Dune includes
#include <config.h>

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>
// #include <dune/quadrature/fixedorder.hh>

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
  // Polynomial and ODE order
  // Grid:
  int N=40;                
  if (argc>1) 
    N = atoi(argv[1]);
  SStruct s(N);
  // GridType grid(s.n_, s.l_, s.h_);
  GridType grid("quadrat.git");
  grid.globalRefine(N);
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
  ModelType model(grid,problem);
  // *** Fluxes 
  FluxType eulerflux(model);
  // *** Operator typedefs
  DgType dg(grid,eulerflux);
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
  ode.printGrid(1, U);
  //typedef DofManager<GridType> DofManagerType;
  //typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;
  //DofManagerType& dm = DofManagerFactoryType :: getDofManager( grid );  
  {
    GrapeDataIO<GridType> dataio;
    std::string gridfile("grid");
    dataio.writeGrid(grid, xdr, gridfile.c_str(), t, n);
    std::string datafile("df");
    dataio.writeData(U, xdr, datafile.c_str(), n);
    //td::string dmfile("dm");
    //dm.write(xdr, dmfile.c_str(),n);
    nextsave+=savestep;
    ++n;
  }
  while (t<1) {
    t=ode.solve(U);
    cout << t << endl;
    ode.printGrid(1, U);
    if (t>nextsave) {
      GrapeDataIO<GridType> dataio;
      std::string gridfile("grid");
      dataio.writeGrid(grid, xdr, gridfile.c_str(), t, n);
      std::string datafile("df");
      dataio.writeData(U, xdr, datafile.c_str(), n);
      //std::string dmfile("dm");
      //dm.write(xdr, dmfile.c_str(),n);
      nextsave+=savestep;
      ++n;
    }
  }
  ode.printGrid(1, U);
} 




