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
#include <sstream>
#include <dune/grid/io/visual/grapedatadisplay.hh>
#include <dune/fem/visual/grape/grapetuple.hh>
#include <dune/common/timer.hh>


using namespace Dune;
using namespace std;

#include "models.hh"
#include "stuff.cc"
#include "robstuff.cc"

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

//#include "grapetuple.hh"


template <class DiscModel,
	  class ODE,class Indicator,class Adapt,class DestinationType>
typename DestinationType::RangeType
solve(DiscModel& model,
      ODE& ode,Indicator& ResiduumErr,
      Adapt* adapt,
      int counter,
      DestinationType& U,DestinationType& V,
      double& t,double& dt) {
  typedef typename DestinationType::RangeType RangeType;
  double t0 = t;
  bool done = true;
  RangeType error;
  do {
    ode.setTime(t0);
    t=ode.solve(U,V,ResiduumErr);
    dt = t - t0;
    error = ResiduumErr.calc(model,adapt,ode,t0,dt);
    //! mark elements and adapt grid
    if (adapt) {
      adapt->param().setTime(t0);
      adapt->param().setTimeStepSize(dt);
      adapt->param().setTimeStepNumber(counter);
      adapt->markRefineEntities();
      adapt->adapt();
      done = (error.one_norm() < adapt->getLocalInTimeTolerance());
      std::cerr << int(done) << " " 
		<< error.one_norm() << " "
		<< adapt->getLocalInTimeTolerance() << std::endl;
    }
  } while (!done);
  ResiduumErr.reset();
  if (adapt) {
    adapt->markCoarsenEntities();
    adapt->adapt();
  }
  return error;
}
template <class DiscModel,
	  class ODE,class Indicator,class Adapt,class DestinationType>
void init(DiscModel& model,
	  ODE& ode,Indicator& ResiduumErr,
	  Adapt& adapt,
	  DestinationType& U,DestinationType& V) {
  typedef typename DestinationType::RangeType RangeType;
  double dt;
  bool done = true;
  RangeType error(0); 
  do {
    U.set(0);
    initialize(model.model().problem(),U);
    dt=ode.solve(U,V,ResiduumErr);
    error =  ResiduumErr.calc(model,adapt,ode,0,dt);
    cout << "START INDICATOR: " << dt << " " 
	 <<  error << " " << flush;
    ode.setTime(0.);
    ResiduumErr.reset();
    if (adapt) {
      adapt->param().setTime(0.0);
      adapt->param().setTimeStepSize(dt);
      adapt->param().setTimeStepNumber(0);
      adapt->markInitRefineEntities();
      cout << " ADAPT " << adapt->getLocalInTimeTolerance() << flush;
      adapt->adapt();
      done = (error.one_norm() < adapt->getLocalInTimeTolerance());
    }
    cout << endl;
  } while (!done);
  U.set(0);
  initialize(model.model().problem(),U);
}

int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  if (argc<2) {
    cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
    exit(EXIT_FAILURE);  
  }
	
  MPISTART
						
  // Polynomial and ODE order
  // Grid:
  GridPtr<GridType> grid(argv[1],MPI_COMM_WORLD);
  GrapeDataIO<GridType> dataio;
	
  int repeats = 1;
  if (argc>2)
    repeats=atoi(argv[2]);	
  int startlevel = 0;
  if (argc>3)
    startlevel=atoi(argv[3]);
  double epsilon = 0.01;
  if (argc>4)
    epsilon=atof(argv[4]);
  int probflag = 0;
  if (argc>5) 
    probflag = atoi(argv[5]);
  int graped = 0;
  if (argc>6)
    graped=atoi(argv[6]);
  bool Hadapt = false;
  double hadapt_tol = 0.;
  if (argc>7) {
    hadapt_tol=atof(argv[7]);
    if (hadapt_tol>1e-20) 
      Hadapt = true;
  }
  bool Padapt = false;
  if (argc>8)
    Padapt=atoi(argv[8]);
  
  // CFL:
  double cfl;
  switch (order) {
  case 0: cfl=0.9;  break;
  case 1: cfl=0.3; break;
  case 2: cfl=0.25;  break;
  case 3: cfl=0.05;  break;
  case 4: cfl=0.09; break;
  }	
  if (argc>9)
    cfl=atof(argv[9]);
	
  //cfl = cfl/3.0;
  //cfl /= 2.0;

  InitialDataType problem(epsilon,probflag,true);
  ModelType model(*grid,problem);
  DiscModelType eulerflux(model);

  string myoutput = "eoc.tex";
  EocOutput eocoutput(myoutput);

  Timer timer;
  int maxit=0;
  double zeit=0.;
  double prevzeit=0.;
  double fehler=0.;
  double prevfehler=0.;
  L1Error<DgType::DestinationType> L1err;
  L1L1Error<ODEType> L1L1err;
  DgType::DestinationType::RangeType err,timeerr,reserr,prevreserr;

  int n=0;
  double nextsave=0.;
  double savestep=0.05;
  double maxtime = problem.endtime();
  int level=0;
	
  for(level=0; level < startlevel ; ++level)
    grid->globalRefine(refStepsForHalf);

  for(int eocloop=0;eocloop < repeats; ++eocloop) {
    // *** Operator typedefs
    DgType dg(*grid,eulerflux);
    ODEType ode(dg,rksteps,cfl); 
    // *** Initial data
    DgType::DestinationType U("U", dg.space());
    DgType::DestinationType V("Unew",dg.space());
    IndType ResiduumErr(dg.part(),order,Padapt);
    OutputType output(&U,ResiduumErr.output());
    TimeDiscParamType * timeDiscParam_ = new TimeDiscParamType (0.0,0.0,0);
    AdaptationType *adaptation_ = 0;
    if (repeats == 1 && Hadapt) {
      adaptation_ = new AdaptationType( dg.part() , *timeDiscParam_, hadapt_tol);
      adaptation_->addAdaptiveFunction(&U,&V,
				       &(ResiduumErr.RT_),
				       &(ResiduumErr.RS_),
      				       &(ResiduumErr.rho_),
				       &(ResiduumErr.lambda_),
				       &(ResiduumErr.maxPol_),
				       &(ResiduumErr.maxPolNew_));
    }
    if (eocloop==0) 
      eocoutput.printInput(problem,*grid,ode,argv[1]);
    
    init(eulerflux,ode,ResiduumErr,adaptation_,U,V);

    GrapeTuple<OutputType>::output(dataio,*grid,0.0,0,"sol",".",output);
    if (graped) {
      GrapeDataDisplay< GridType > grape(dg.part());
      grape.addData(ResiduumErr.RT_,"El-Res",-4);
      grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
      grape.addData(ResiduumErr.rho_,"rho",-2);
      grape.addData(ResiduumErr.lambda_,"lam",-1);
      grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
      grape.dataDisplay(U);
    }

    double t=0.0;
    int counter=0;
    int outputcounter=0;
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
      double ldt;
      reserr += solve(eulerflux,ode,ResiduumErr,adaptation_,counter,
		      U,V,t,ldt);
      U.assign(V);
      ++counter;
      if (repeats==1 && counter%50==0) {
	outputcounter++;
        GrapeTuple<OutputType>::
	  output(dataio,*grid,t,outputcounter,"sol",".",output);
        if (graped) {
	  GrapeDataDisplay< GridType > grape(dg.part());
	  grape.addData(ResiduumErr.RT_,"El-Res",-4);
	  grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
	  grape.addData(ResiduumErr.rho_,"rho",-2);
	  grape.addData(ResiduumErr.lambda_,"lam",-1);
	  grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
	  grape.dataDisplay(U);
	}
      }
      timeerr += L1L1err.norm(problem,ode,t-ldt,ldt);
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      if (repeats==1) 
	std::cout << counter << " " << t << " " << ldt << std::endl;
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
    eocoutput.printTexAddError(fehler,prevfehler,reserr.one_norm(),prevreserr.one_norm(),zeit,grid->size(0),counter,averagedt);
    
    outputcounter++;
    GrapeTuple<OutputType>::output(dataio,*grid,t,outputcounter,"sol",".",output);
    if(graped) {
      GrapeDataDisplay< GridType > grape(dg.part());
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
    prevreserr = reserr;
    ++maxit;
  }
  
  eocoutput.printTexEnd(timer.elapsed());
  
  MPIEND
} 




