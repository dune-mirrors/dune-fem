bool Hadapt = false;
bool Padapt = false;
// Dune includes
#include <config.h>

#include <dune/fem/pass/utility.hh>

#include <dune/common/misc.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/space/dgspace/dgleafindexset.hh>

#include <dune/fem/io/file/grapedataio.hh>
#include <dune/fem/space/common/boundary.hh>

#include <iostream>
#include <string>

#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/common/timer.hh>

using namespace Dune;
using namespace std;

#include "models.hh"
#include "stuff.cc"
#include "robstuff.cc"
#include "aposteriori.hh"

// include adaptation interface class
#include "adaptation.hh"

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

int main(int argc, char ** argv, char ** envp) {

  // *** Typedefs for adaptation 
  typedef TimeDiscrParam                           TimeDiscParamType;

  typedef Adaptation <DgType::DestinationType, 
          TimeDiscParamType>                       AdaptationType;

  typedef DgType::DestinationType::DiscreteFunctionSpaceType::GridPartType GridPartType;
  typedef GridPartType::IndexSetType IndexSetType;

  // *** Initialization
  if (argc<2) {
    cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
    exit(EXIT_FAILURE);  
  }
	
  MPISTART
						
  // Polynomial and ODE order
  // Grid:
  GridPtr<GridType> grid(argv[1],MPI_COMM_WORLD);
	
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
	

  bool Hadapt = false;
  if (argc>6)
    Hadapt=atoi(argv[6]);
  bool Padapt = false;
  if (argc>7)
    Padapt=atoi(argv[7]);
  
  // CFL:
  double cfl;
  switch (order) {
  case 0: cfl=0.9;  break;
  case 1: cfl=0.2; break;
  case 2: cfl=0.1;  break;
  case 3: cfl=0.05;  break;
  case 4: cfl=0.09; break;
  }
	
  if (argc>8)
    cfl=atof(argv[8]);
	
  cfl = cfl/3.0;
  cfl /= 2.0;

  cout << epsilon << endl;
	
  InitialDataType problem(epsilon,true);
	
  string myoutput = "eoc.tex";
  EocOutput eocoutput(myoutput);
  // *** Models
  ModelType model(*grid,problem);
  // *** Fluxes 
  DiscModelType eulerflux(model);

  Timer timer;
  int maxit=0;
  double zeit=0.;
  double prevzeit=0.;
  double fehler=0.;
  double prevfehler=0.;
	
  L1Error<DgType::DestinationType> L1err;
  L1L1Error<ODEType> L1L1err;
  DgType::DestinationType::RangeType err,timeerr,reserr;

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
    DgType dg(*grid,eulerflux); // ,upwind);
    ODEType ode(dg,rksteps,cfl); 
    // *** Initial data
    DgType::DestinationType U("U", dg.space());
    DgType::DestinationType tmp("tmp",dg.space());
    Residuum<GridType,DiscModelType,ODEType> ResiduumErr(*grid,Padapt);
    IndexSetType * iset = new IndexSetType ( *grid );
    GridPartType * gridPart_ = new GridPartType ( *grid, *iset );
    // initialize time discretization parameters
    TimeDiscParamType * timeDiscParam_ = new TimeDiscParamType (0.0,0.0,0);
    // initialize adaptation if wanted
    const char * paramfile = 0;
    AdaptationType *adaptation_ = 0;
    if (repeats == 1 && Hadapt) {
      AdaptationType *adaptation_ = 
	new AdaptationType( *gridPart_ , *timeDiscParam_, paramfile);
      adaptation_->addAdaptiveFunction(&U);
    }

    if (eocloop==0) 
      eocoutput.printInput(problem,*grid,ode,argv[1]);
    for (int i=0;i<5;++i) {
      U.set(0);
      initialize(problem,U);
      // dg.limit(U,tmp);
      double dt=ode.solve(U,ResiduumErr);
      cout << "START INDICATOR: " << dt << " " 
	   << ResiduumErr.calc(eulerflux,adaptation_,ode,0,dt) << endl;
      ode.setTime(0.);
    }
    if (graped) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.addData(ResiduumErr.RT_,"El-Res",-4);
      grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
      grape.addData(ResiduumErr.rho_,"rho",-2);
      grape.addData(ResiduumErr.lambda_,"lam",-1);
      grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
      grape.dataDisplay(U);
    }
    U.set(0);
    initialize(problem,U);
    dg.limit(U,tmp);
    ode.setTime(0.);
    double t=0.0;
    int counter=0;
    FieldVector<double,1> projectionError = L1err.norm(problem,U,t);  
    cout << "Projection error " << 
      problem.myName << ": " << projectionError << endl;
	
    double maxdt=0.,mindt=1.e10,averagedt=0.;
    timeerr = 0.;
    err = 0.;
    reserr = 0.;
    
    // *** Time loop
    while (t<maxtime) 
    {
      double ldt = -t;
      t=ode.solve(U,ResiduumErr);
      ldt += t;
      reserr  += ResiduumErr.calc(eulerflux,adaptation_,
				  ode,t-ldt,ldt);
      timeDiscParam_->setTime(t);
      timeDiscParam_->setTimeStepSize(ldt);
      timeDiscParam_->setTimeStepNumber(counter);
      
      //! mark elements and adapt grid
      if (adaptation_) {
	adaptation_->markEntities();
	adaptation_->adapt();
      }
      //initialize(problem,U);
      if (graped && repeats==1 && Hadapt) {
	GrapeDataDisplay< GridType > grape(*grid);
	grape.dataDisplay(U);
      }
     
      timeerr += L1L1err.norm(problem,ode,t-ldt,ldt);
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      dg.limit(U,tmp);
      dg.switchupwind();
      ++counter;
    }
    
    averagedt /= double(counter);

    err = L1err.norm(problem,U,t);
    cout << "Error " << problem.myName << ": " << err << endl;
    cout << "TimeSpace-Error " << 
      problem.myName << ": " << timeerr << endl;
    cout << "Residuum-Error " << 
      problem.myName << ": " << reserr << endl;
    cout << endl;
    fehler = err;		
    zeit = timer.elapsed()-prevzeit;
    eocoutput.printTexAddError(fehler,prevfehler,zeit,grid->size(0),counter,averagedt);
    
    if(graped){ // && eocloop == repeats-1) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.addData(ResiduumErr.RT_,"El-Res",-4);
      grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
      grape.addData(ResiduumErr.rho_,"rho",-2);
      grape.addData(ResiduumErr.lambda_,"lam",-1);
      grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
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




