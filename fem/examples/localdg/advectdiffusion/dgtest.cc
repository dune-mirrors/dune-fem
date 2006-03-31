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

#include <dune/io/visual/grapedatadisplay.hh>
#include <dune/common/timer.hh>

using namespace Dune;
using namespace std;

#include "models.hh"
#include "stuff.cc"
#include "robstuff.cc"

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  if (argc<2) {
    cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
    exit(EXIT_FAILURE);  
  }
	
  MPISTART
						
  // Polynomial and ODE order
  // Grid:
  GridType* grid=MacroGrid<GridType>().generate(argv[1],MPI_COMM_WORLD);
	
  int repeats = 1;
  if (argc>2)
    repeats=atoi(argv[2]);
	
  int startlevel = 0;
  if (argc>3)
    startlevel=atoi(argv[3]);
	
  double epsilon = 0.01;
  if (argc>4)
    epsilon=atof(argv[4]);
	
  int graped = 0;
  if (argc>5)
    graped=atoi(argv[5]);
		
  // CFL:
  double cfl;
  switch (order) {
  case 0: cfl=0.9;  break;
  case 1: cfl=0.2; break;
  case 2: cfl=0.1;  break;
  case 3: cfl=0.05;  break;
  case 4: cfl=0.09; break;
  }
	
  if (argc>6)
    cfl=atof(argv[6]);
	
  cfl = cfl/3.0;
  cfl /= 2.0;

  cout << epsilon << endl;
	
  InitialDataType problem(epsilon,true);
	
  string myoutput = "eoc.tex";
  EocOutput eocoutput(myoutput);
  // *** Models
  ModelType model(*grid,problem);
  // *** Fluxes 
  FluxType eulerflux(model);

  Timer timer;
  int maxit=0;
  double zeit=0.;
  double prevzeit=0.;
  double fehler=0.;
  double prevfehler=0.;
	
  L2Error<DgType::DestinationType> L2err;
  FieldVector<double,1> err;

  // printSGrid(0,0,dg.space(),U);
	
  int n=0;
  double nextsave=0.;
  double savestep=0.05;
  double maxtime = problem.endtime();
  int level=0;
	
  FieldVector<double,dimworld> upwind(1.);
  upwind[0] *= 0.37;
	
  for(level=0; level < startlevel ; ++level)
    grid->globalRefine(refStepsForHalf);

  for(int eocloop=0;eocloop < repeats; ++eocloop) {
    // *** Operator typedefs
    DgType dg(*grid,eulerflux,upwind);
    ODEType ode(dg,rksteps,cfl); 
    // *** Initial data
    DgType::DestinationType U("U", dg.space());
    DgType::DestinationType tmp("tmp",dg.space());
    // DofManagerType& dm = DofManagerFactoryType :: getDofManager( *grid );  
    // grid->preAdapt();
    // dm.resizeForRestrict();
    // grid->globalRefine(refStepsForHalf);
    // dm.resize();
    // dm.dofCompress();
    // grid->postAdapt();
    
    /*		 DgType::DestinationType tmp("tmp",dg.space());
      {
      initialize(*problem,U);
      dg.limit(U,tmp);
      }
      ode.solve(U);
      dg.limit(U,tmp);
      
    */
    
    initialize(problem,U);
    dg.limit(U,tmp);

    if (eocloop==0) 
      eocoutput.printInput(problem,grid,ode,argv[1]);
    
    if(graped) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.dataDisplay(U);
    }
    
    double t=0.0;
    int counter=0;
    FieldVector<double,1> projectionError = L2err.norm(problem,U,t);  
    cout << "Projection error " << problem.myName << ": " << projectionError << endl;
	
    double maxdt=0.,mindt=1.e10,averagedt=0.;
    // *** Time loop
    while (t<maxtime) {
      double ldt = -t;
      t=ode.solve(U);
      ldt += t;
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      dg.limit(U,tmp);
      dg.switchupwind();
      if(counter%100 == 0) {
	err = L2err.norm(problem,U,t);
	if(err > 1e5 || ldt < 1e-10) {
	  averagedt /= double(counter);
	  cout << "Solution doing nasty things!" << std::endl;
	  cout << t << endl;
	  /*{
	    GrapeDataDisplay< GridType > grape(*grid);
	    grape.dataDisplay(U);
	    }*/
	  eocoutput.printTexAddError(err[0],prevfehler,-1,grid->size(0),counter,averagedt);
	  eocoutput.printTexEnd(timer.elapsed());
	  exit(EXIT_FAILURE);
	}
      }
      ++counter;
    }
    
    averagedt /= double(counter);

    err = L2err.norm(problem,U,t);
    cout << "Error " << problem.myName << ": " << err << endl;
    
    fehler = err;		
    zeit = timer.elapsed()-prevzeit;
    eocoutput.printTexAddError(fehler,prevfehler,zeit,grid->size(0),counter,averagedt);
    
    if(graped && eocloop == repeats-1) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.dataDisplay(U);
    }
    
    if(zeit > 3000.)
      break;
    
    if(eocloop < repeats-1) {
      grid->globalRefine(refStepsForHalf);
      ++level;
    }
    prevzeit = zeit;
    prevfehler = fehler;
    ++maxit;
  }
  
  eocoutput.printTexEnd(timer.elapsed());
  
  MPIEND
} 




