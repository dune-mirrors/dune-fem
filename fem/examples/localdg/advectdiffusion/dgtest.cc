// Dune includes
#include <config.h>

#if PROBLEM == 5 
#define EULER_PERFORMANCE
const double globalTimeStep = 1e-3;
#endif

// inlcude template meta fvector 
#include "fvector.hh"

#include <dune/common/utility.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/common/misc.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/io/file/grapedataio.hh>

#include <iostream>
#include <string>

#if HAVE_GRAPE
#include <dune/grid/io/visual/grapedatadisplay.hh>
#endif

#include <dune/common/timer.hh>

#include <dune/common/mpihelper.hh>

using namespace Dune;
using namespace std;

#include "models.hh"
#include "stuff.cc"
#include "robstuff.cc"

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

int main(int argc, char ** argv, char ** envp) {

  MPIHelper::instance(argc,argv);

  try {

  // *** Initialization
  if (argc<2) {
    cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
    exit(EXIT_FAILURE);  
  }
	
  // Polynomial and ODE order
  // Grid:
  GridPtr<GridType> grid(argv[1]); // ,MPI_COMM_WORLD);
	
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
  switch (order) 
  {
    case 0: cfl=0.9;  break;
    case 1: cfl=0.2; break;
    case 2: cfl=0.15;  break;
    case 3: cfl=0.05;  break;
    case 4: cfl=0.09; break;
  }

//#if PROBLEM == 5 
  if (argc>6)
    cfl=atof(argv[6]);
//#endif
  
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
  FieldVector<double,ModelType::dimRange> err;

  // printSGrid(0,0,dg.space(),U);
	
  //int n=0;
  //double nextsave=0.;
  //double savestep=0.05;
  double maxtime = problem.endtime();
  int level=0;
	
  FieldVector<double,dimworld> upwind(1.);
  upwind[0] *= 0.37;
	
  for(level=0; level < startlevel ; ++level)
    grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());

  for(int eocloop=0;eocloop < repeats; ++eocloop) {
    // *** Operator typedefs
    DgType dg(*grid,eulerflux,upwind);
    ODEType ode(dg,rksteps,cfl); 

    // set timestep size as it is fixed 
#ifdef EULER_PERFORMANCE
    ode.provideTimeStepEstimate(globalTimeStep);
    ode.syncTimeStep();
#endif
    
    // *** Initial data
    DgType::DestinationType U("U", dg.space());
    DgType::DestinationType tmp("tmp",dg.space());
    // DofManagerType& dm = DofManagerFactoryType :: getDofManager( *grid );  
    // grid->preAdapt();
    // dm.resizeForRestrict();
    // grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
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
      eocoutput.printInput(problem,*grid,ode,argv[1]);
    
#if HAVE_GRAPE 
    if(graped > 0) {
      GrapeDataDisplay< GridType > grape(U.space().gridPart());
      grape.dataDisplay(U);
    }
#endif 
    
    double t=0.0;
    int counter=0;
    FieldVector<double,ModelType::dimRange> projectionError = L2err.norm(problem,U,t);  
    cout << "Projection error " << problem.myName << ": " << projectionError << endl;
	
    double maxdt=0.,mindt=1.e10,averagedt=0.;
    // *** Time loop
    while (t<maxtime) 
    {
      double ldt = -t;
      t=ode.solve(U);

#ifdef EULER_PERFORMANCE
      // only for Euler Performance Test
      ode.resetTimeStepEstimate();
      ode.provideTimeStepEstimate(globalTimeStep);
      ode.syncTimeStep();
#endif
      
#if HAVE_GRAPE 
      if(graped > 0)
      {
        if(counter%graped == 0 && counter > 0) 
        {
          GrapeDataDisplay< GridType > grape(U.space().gridPart());
          grape.dataDisplay(U);
        }
      }
#endif 
      ldt += t;
      //std::cout << "Current time = " << t << "\n";
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      // dg.limit(U,tmp);
      dg.switchupwind();
      if(0 && counter%100 == 0) 
      {
        err = L2err.norm(problem,U,t);
        if(err.one_norm() > 1e5 || ldt < 1e-10) 
        {
          averagedt /= double(counter);
          cout << "Solution doing nasty things!" << std::endl;
          cout << t << endl;
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
    
    fehler = err.two_norm();		
    zeit = timer.elapsed()-prevzeit;
    eocoutput.printTexAddError(fehler,prevfehler,zeit,grid->size(0),counter,averagedt);
   
#if HAVE_GRAPE
    if(graped>0 && eocloop == repeats-1) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.dataDisplay(U);
    }
#endif
    
    if(zeit > 3000.)
      break;
    
    if(eocloop < repeats-1) {
      grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());
      ++level;
    }
    prevzeit = zeit;
    prevfehler = fehler;
    ++maxit;
  }
  
  eocoutput.printTexEnd(timer.elapsed());

  }
  catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "Generic exception!" << std::endl;
    return 2;
  }

  return 0;  
} 




