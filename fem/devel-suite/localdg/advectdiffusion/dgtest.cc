// Dune includes
#include <config.h>

#if PROBLEM == 5 
#define EULER_PERFORMANCE
const double globalTimeStep = 1e-3;
#endif

#include <dune/common/fvector.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/grid/common/gridpart.hh>

#include <dune/common/misc.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/fem/misc/l2error.hh>

#include <iostream>
#include <string>

#include <dune/common/timer.hh>
#include <dune/common/mpihelper.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/io/file/datawriter.hh>

#include "models.hh"
#include "stuff.cc"

using namespace std;
using namespace Dune;

typedef DofManager<GridType> DofManagerType;
typedef DofManagerFactory<DofManagerType> DofManagerFactoryType;

int main(int argc, char ** argv, char ** envp) {

  MPIHelper::instance(argc,argv);

  try {

  // *** Initialization
	
  Parameter::append(argc, argv);
  if (argc==2) Parameter::append(argv[1]);
  else Parameter::append("parameter");
    
  // Polynomial and ODE order
  // Grid:
	
  std::string filename;
  Parameter::get("fem.localdg.gridfile", filename);
  GridPtr<GridType> grid(filename); // ,MPI_COMM_WORLD);
  double startTime = Parameter::getValue<double>("fem.localdg.starttime",0.0);
  int repeats = Parameter::getValue<int>("fem.localdg.repeats");
  int startlevel = Parameter::getValue<int>("fem.localdg.startlevel");
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
  Parameter::get("fem.localdg.cfl", cfl, cfl);
  std::cout << " CFL : " << cfl << std::endl;

  int printCount = Parameter::getValue("fem.localdg.printcount",-1);

  InitialDataType problem;
	
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

  double endTime = problem.endtime();
  int level=0;
	
  FieldVector<double,dimworld> upwind(1.);
  upwind[0] *= 0.37;
	
  for(level=0; level < startlevel ; ++level)
    grid->globalRefine(DGFGridInfo<GridType>::refineStepsForHalf());

  for(int eocloop=0;eocloop < repeats; ++eocloop) {
    // *** Operator typedefs
    DgType dg(*grid,eulerflux,upwind);
    TimeProvider sertp(startTime,cfl);
    typedef GridType :: Traits :: CollectiveCommunication CommunicatorType;
    ParallelTimeProvider<CommunicatorType> tp(grid->comm(),sertp);
    ODEType ode(dg,sertp,rksteps,Parameter::verbose()); 
    dg.timeProvider(&sertp);
    
    // *** Initial data
    DgType::DestinationType U("U", dg.space());
    DgType::DestinationType tmp("tmp",dg.space());
    
    initialize(problem,U);
    dg.limit(U,tmp);

    if (eocloop==0) 
      eocoutput.printInput(problem,*grid,ode,argv[1]);
    
    IOTupleType dataTup ( &U );
    typedef DataWriter< GridType, IOTupleType > DataWriterType;
    DataWriterType dataWriter( *grid, filename, dataTup,startTime,endTime);
    
    double t=startTime;
    int counter=0;
    FieldVector<double,ModelType::dimRange> projectionError = L2err.norm(problem,U,t);  
    cout << "Projection error " << problem.myName << ": " << projectionError << endl;


    const double globalTimeStep = Parameter::getValue<double>("fem.localdg.global_tstep");
		double maxdt=0.,mindt=1.e10,averagedt=0.;
    // *** Time loop
    dataWriter.write(t , counter );
    tp.resetTimeStepEstimate(globalTimeStep);
    ode.initialize(U);
    tp.syncTimeStep();

    while (t<endTime) 
    {
      tp.resetTimeStepEstimate(globalTimeStep);
      ode.solve(U);
      tp.augmentTime();
      tp.syncTimeStep();
      t = tp.time();
      double ldt = tp.deltaT();
      if (!U.dofsValid()) {
	      std::cout << "Invalid DOFs" << std::endl;
	      dataWriter.write(1e10, counter );  
	      abort();
      }
      
      if (printCount>0 && counter%printCount==0) {
        std::cout << "step: " << counter << " time: " << t << " deltaT:" << ldt << std::endl;
      }
      dataWriter.write(t, counter );  
      
      mindt = (ldt<mindt)?ldt:mindt;
      maxdt = (ldt>maxdt)?ldt:maxdt;
      averagedt += ldt;
      dg.limit(U,tmp);
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
    
    dataWriter.write(t, counter );  
    averagedt /= double(counter);

    err = L2err.norm(problem,U,t);
    cout << "Error " << problem.myName << ": " << err << endl;
    
    fehler = err.two_norm();		
    zeit = timer.elapsed()-prevzeit;
    eocoutput.printTexAddError(fehler,prevfehler,zeit,grid->size(0),counter,averagedt);
   
    
    
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
  Parameter::write("parameter.log");

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




