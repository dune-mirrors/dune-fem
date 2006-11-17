#undef HAVE_MPI
#undef HAVE_MPI_CPP
#include <dune/grid/io/file/dgfparser/gridtype.hh>
#include <iostream>
#include <string>
#include <sstream>
#include <vector> 
// #include <dune/grid/yaspgrid.hh>
// #include <dune/grid/sgrid.hh>
// #include <dune/grid/albertagrid.hh>
#include "advectdiff.hh"

#include <dune/fem/solver/odesolver.hh>

// Approximations Ordnung
enum {order=POLORDER,rksteps=POLORDER+1}; 
// enum {order=POLORDER,rksteps=POLORDER}; 
// Gitterauswahl

// Modell- und Flussauswahl
// Skalar
#if PROBLEM == 1
   #include "scalarmodels.hh"
   #include "initadvectdiff.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> DiscModelType;
   typedef UpwindFlux<ModelType> DiscModelType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 2
   #include "scalarmodels.hh"
   #include "initburgers.cc"
   typedef U0<GridType> InitialDataType;
   typedef BurgersModel<GridType,InitialDataType > ModelType;
   typedef LLFFlux<ModelType> DiscModelType;
   //typedef DGLimitedAdvectionOperator<ModelType,LLFFlux,order> DgType;
   typedef DGAdvectionOperator<ModelType,LLFFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 3
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> DiscModelType;
   typedef UpwindFlux<ModelType> DiscModelType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 4
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> DiscModelType;
   typedef UpwindFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 5
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> DiscModelType;
   typedef UpwindFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 6
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0RotCone<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> DiscModelType;
   typedef UpwindFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 7
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0BuckLev<GridType> InitialDataType; 
   typedef BuckLevModel<GridType,InitialDataType > ModelType;
   typedef LFFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,LFFlux,order> DgType;
   // typedef UpwindFlux<ModelType> DiscModelType;
   // typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 11
#include "euler_mhd/eulermodel.hh"
   typedef U0Sod InitialDataType; 
   typedef EulerModel<GridType,InitialDataType > ModelType;
   //typedef DWNumFlux<ModelType> DiscModelType;
   typedef LFFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,LFFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#elif PROBLEM == 12
#include <dune/fem/examples/localdg/advectdiffusion/euler_mhd/eulermodel.hh>
   typedef FFS InitialDataType; 
   typedef EulerModel<GridType,InitialDataType > ModelType;
   //typedef DWNumFlux<ModelType> DiscModelType;
   typedef LFFlux<ModelType> DiscModelType;
   typedef DGAdvectionOperator<ModelType,LFFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType,rksteps> ODEType;
#endif
#include "aposteriori.hh"
#include "adaptation.hh"
typedef Residuum<DgType::GridPartType,DiscModelType,ODEType> IndType; 
typedef Pair<DgType::DestinationType*,IndType::OutputType> OutputType;
typedef TimeDiscrParam TimeDiscParamType;
typedef Adaptation <DgType::DestinationType,TimeDiscParamType>  AdaptationType;
