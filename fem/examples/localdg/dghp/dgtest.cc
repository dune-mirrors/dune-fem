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
#include <sstream>
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

template <class T1,class T2,int N>
struct GrapeTupleHelper {
  template <class DataIO>
  static void apply(DataIO& dataio,string name,int n,
	     Pair<T1,T2>& tup) {
    stringstream dataname;
    dataname << name << "_" << N;
    dataio.writeData(*(tup.first()), xdr, dataname.str().c_str(), n);
    GrapeTupleHelper<typename T2::Type1,typename T2::Type2,N+1>::apply(dataio,name,n,tup.second());
  }
};
template <class T1,int N>
struct GrapeTupleHelper<T1,Nil,N> {
  template <class DataIO>
  static void apply(DataIO& dataio,string name,int n,
      Pair<T1,Nil>& tup) {
    stringstream dataname;
    dataname << name << "_" << N;
    dataio.writeData(*(tup.first()), xdr, dataname.str().c_str(), n);    
  }
};
struct GrapeTupleOutput {
  template <class DataIO,class GridType,class Tup>
  static void apply(DataIO& dataio,GridType& grid,double t,int n,
	     string name,Tup& tup) {
    string gridname = name + "_grid";
    dataio.writeGrid(grid, xdr, gridname.c_str(), t, n);
      GrapeTupleHelper<typename Tup::Type1,typename Tup::Type2,0>::
      apply(dataio,name,n,tup);
  }
};

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
    t=ode.solve(U,V,ResiduumErr);
    dt = t - t0;
    error = ResiduumErr.calc(model,adapt,ode,t0,dt);
    //! mark elements and adapt grid
    if (adapt) {
      adapt->param().setTime(t0);
      adapt->param().setTimeStepSize(dt);
      adapt->param().setTimeStepNumber(counter);
      adapt->markEntities();
      adapt->refine();
      done = (error.one_norm() < adapt->getLocalTolerance());
    }
  } while (!done);
  return error;
}
template <class DiscModel,
	  class ODE,class Indicator,class Adapt,class DestinationType>
void init(DiscModel& model,
	  ODE& ode,Indicator& ResiduumErr,
	  Adapt& adapt,
	  DestinationType& U,DestinationType& V) {
  for (int i=0;i<5;++i) {
    U.set(0);
    initialize(model.model().problem(),U);
    double dt=ode.solve(U,V,ResiduumErr);
    cout << "START INDICATOR: " << dt << " " 
	 << ResiduumErr.calc(model,adapt,ode,0,dt) << endl;
    ode.setTime(0.);
  }
  U.set(0);
  initialize(model.model().problem(),U);
}

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
	
  //cfl = cfl/3.0;
  //cfl /= 2.0;

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
    DgType::DestinationType V("Unew",dg.space());
    typedef Residuum<GridType,DiscModelType,ODEType> IndType; 
    IndType ResiduumErr(*grid,Padapt);
    Pair<DgType::DestinationType*,IndType::OutputType> 
      output(&U,ResiduumErr.output());
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
      adaptation_->addAdaptiveFunction(&V);
      adaptation_->addAdaptiveFunction(&(ResiduumErr.RT_));
      adaptation_->addAdaptiveFunction(&(ResiduumErr.RS_));
      adaptation_->addAdaptiveFunction(&(ResiduumErr.rho_));
      adaptation_->addAdaptiveFunction(&(ResiduumErr.lambda_));
      adaptation_->addAdaptiveFunction(&(ResiduumErr.maxPol_));
    }

    if (eocloop==0) 
      eocoutput.printInput(problem,*grid,ode,argv[1]);
    
    init(eulerflux,ode,ResiduumErr,adaptation_,U,V);

    if (graped) {
      GrapeDataDisplay< GridType > grape(*grid);
      grape.addData(ResiduumErr.RT_,"El-Res",-4);
      grape.addData(ResiduumErr.RS_,"Jmp-Res",-3);
      grape.addData(ResiduumErr.rho_,"rho",-2);
      grape.addData(ResiduumErr.lambda_,"lam",-1);
      grape.addData(ResiduumErr.maxPol_,"poldeg",-1);
      grape.dataDisplay(U);
    }

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
      double ldt;
      reserr += solve(eulerflux,ode,ResiduumErr,adaptation_,counter,
		      U,V,t,ldt);
      U.assign(V);
      //initialize(problem,U);
      ++counter;
      if (repeats==1 && counter%1000==0) {
	//GrapeDataDisplay< GridType > grape(*grid);
	//grape.dataDisplay(U);
	GrapeTupleOutput::apply(dataio,*grid,t,counter/200,"output",output);
      }
     
      timeerr += L1L1err.norm(problem,ode,t-ldt,ldt);
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
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




