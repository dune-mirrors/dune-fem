// Dune includes
#include <config.h>

#include <dune/fem/pass/utility.hh>

#include <dune/common/misc.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/space/dgspace/dgleafindexset.hh>

#include <dune/fem/io/file/grapedataio.hh>
#include <dune/fem/space/common/boundary.hh>

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
      double& t,double& dt,
      double& locerror,
      bool start=false) {
  typedef typename DestinationType::RangeType RangeType;
  double t0 = t;
  bool done = true;
  RangeType error;
  L1L1Error<ODE> L1L1err;
  locerror = 0.;
  do {
    ode.setTime(t0);
    t=ode.solve(U,V,ResiduumErr);
    dt = t - t0;
    // locerror = L1L1err.norm(model.model().problem(),ode,t0,dt);
    error = ResiduumErr.calc(model,adapt,ode,t0,dt);
    //! mark elements and adapt grid
    if (adapt) {
      int old_size = ResiduumErr.part_.grid().size(0);
      adapt->param().setTime(t0);
      adapt->param().setTimeStepSize(dt);
      adapt->param().setTimeStepNumber(counter);
      adapt->markRefineEntities();
      adapt->adapt();
      std::cout << "   " << "loop:" << int(done) << " " << counter<< " "
		<< error.one_norm() << " "
		<< adapt->getLocalInTimeTolerance() << " " 
		<< ResiduumErr.numPAdapt() << "    "
		<< old_size << " " << ResiduumErr.part_.grid().size(0) << " " 
		<< std::endl;
      done = (error.one_norm() < adapt->getLocalInTimeTolerance());
    }
  } while (!done);
  ResiduumErr.reset();
  ode.project(V,ResiduumErr);
  if (!start && adapt) {
    adapt->param().setTime(t0);
    adapt->param().setTimeStepSize(dt);
    adapt->param().setTimeStepNumber(counter);
    adapt->markCoarsenEntities(V,ResiduumErr.rho_);
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
  int padapt_num = 0;
  do {
    U.set(0);
    initialize(model.model().problem(),U);
    dt=ode.solve(U,V,ResiduumErr);
    error =  ResiduumErr.calc(model,adapt,ode,0,dt,true);
    cout << "START INDICATOR: " << dt << " " 
	 <<  error << " " << flush;
    ode.setTime(0.);
    if (adapt) {
      adapt->param().setTime(0);
      adapt->param().setTimeStepSize(dt);
      adapt->param().setTimeStepNumber(0);
      adapt->markInitRefineEntities();
      done = (error.one_norm() < adapt->getInitTolerance());
      std::cout << int(done) << " " 
		<< error.one_norm() << " "
		<< adapt->getInitTolerance() << std::endl;
      adapt->adapt();
    } else {
      done = (padapt_num >= ResiduumErr.numPAdapt());
      padapt_num = ResiduumErr.numPAdapt();
      std::cout << "number of p-adaptive elements " << padapt_num << std::endl;
    }
    cout << endl;
    ResiduumErr.reset();
  } while (!done);
  
  ode.setTime(0.);
  U.set(0);
  initialize(model.model().problem(),U);
  if (adapt) {
    ode.project(U,ResiduumErr);
    for (int i=0;i<15;i++) {
      double t = 0.0,locerror;
      solve(model,ode,ResiduumErr,adapt,0,U,V,t,dt,locerror);
      ode.setTime(0.);
      U.set(0);
      initialize(model.model().problem(),U);
      ode.project(U,ResiduumErr);
    }
  }
}

int main(int argc, char ** argv, char ** envp) {
  // *** Initialization
  if (argc<2) {
    cout << "Call: dgtest gridfilename [ref-steps=1] [start-level=0] [epsilon=0.01] [use-grape=0]" << endl;
    exit(EXIT_FAILURE);  
  }
	
  
  // Polynomial and ODE order
  // Grid:
  GridPtr<GridType> grid(argv[1]);
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
    hadapt_tol=atof(argv[7])*1.;
    if (hadapt_tol>1e-20) 
      Hadapt = true;
  }
  bool Padapt = false;
  if (argc>8)
    Padapt=atoi(argv[8]);
  
  // CFL:
  double cfl;
  switch (order) {
  case 0: cfl=0.475;  break;
  case 1: cfl=0.3; break;
  case 2: cfl=0.24;  break;
  case 3: cfl=0.19;  break;
  case 4: cfl=0.15; break;
  }	
  if (argc>9)
    cfl=atof(argv[9]);

  // cfl *= 2.;
	
  InitialDataType problem(epsilon,probflag,false);
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
  double savestep = problem.saveinterval();
  double nextsave = savestep;
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
      adaptation_ = new AdaptationType( dg.part() , *timeDiscParam_, hadapt_tol,
					problem.endtime());
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
    if (repeats>1) {
      GrapeTuple<OutputType>::output(dataio,*grid,0.0,eocloop*10,
				     "sol",".",output);
    } else {
      GrapeTuple<OutputType>::output(dataio,*grid,0.0,0,"sol",".",output);
    }
    if (graped) {
      GrapeDataDisplay< GridType > grape(dg.part());
      grape.addData(ResiduumErr.RT_,"El-Res",-4);
      grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
      grape.addData(ResiduumErr.rho_,"rho",-2);
      grape.addData(ResiduumErr.lambda_,"lam",-1);
      grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
      grape.dataDisplay(U);
    }

    double t=0.;
    int counter=0;
    int outputcounter=0;
    FieldVector<double,1> projectionError = L1err.norm(problem,U,t).two_norm();  
    cout << "Projection error " << 
      problem.myName << ": " << projectionError << endl;
	
    double maxdt=0.,mindt=1.e10,averagedt=0.;
    timeerr = 0.;
    err = 0.;
    reserr = 0.;

    int totalsize = grid->size(0);
    
    // *** Time loop
    while (t<maxtime) 
    {
      double ldt;
      double locerror;
      reserr += solve(eulerflux,ode,ResiduumErr,adaptation_,counter,
		      U,V,t,ldt,locerror);
      // t=ode.solve(U,V,ResiduumErr);
      U.assign(V);
      totalsize += grid->size(0);
      ++counter;
      if (repeats==1 && t>nextsave) {
        nextsave += savestep;
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
      } else if (repeats>1 && t-ldt<maxtime*0.25 && t>maxtime*0.25) {
	 GrapeTuple<OutputType>::output(dataio,*grid,0.0,eocloop*10+1,
				     "sol",".",output);
      } 
      timeerr += locerror;
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      if (repeats==1) 
	std::cout << counter << " " << t << " " << ldt << " "
		  << locerror/ldt << " " << reserr << " "
		  << ResiduumErr.numPAdapt() << " "
		  << "# time-steps" << std::endl;
    }
    
    averagedt /= double(counter);

    err = L1err.norm(problem,U,t);
    cout << "Error " << problem.myName << ": " << err << endl;
    cout << "TimeSpace-Error " << 
      problem.myName << ": " << timeerr << endl;
    cout << "Residuum-Error " << 
      problem.myName << ": " << reserr << endl;
    cout << endl;
    fehler = err.two_norm();		
    zeit = timer.elapsed()-prevzeit;
    eocoutput.printTexAddError(fehler,prevfehler,
			       reserr.one_norm(),prevreserr.one_norm(),
			       zeit,totalsize,counter,averagedt);
    
    outputcounter++;
    if (repeats>1) {
      GrapeTuple<OutputType>::output(dataio,*grid,0.0,eocloop*10+2,
				     "sol",".",output);
    }     
    else 
      GrapeTuple<OutputType>::
	output(dataio,*grid,t,outputcounter,"sol",".",output);
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
  
} 




