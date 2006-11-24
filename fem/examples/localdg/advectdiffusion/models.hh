#undef HAVE_MPI
#include <dune/grid/io/file/dgfparser/gridtype.hh>
// #include <dune/grid/yaspgrid.hh>
// #include <dune/grid/sgrid.hh>
// #include <dune/grid/albertagrid.hh>
#include "advectdiff.hh"

#include <dune/fem/solver/odesolver.hh>
// #include "odesolver.hh"
// Approximations Ordnung
enum {order=POLORDER,rksteps=POLORDER+1}; 
// Gitterauswahl

// Modell- und Flussauswahl
// Skalar
#if PROBLEM == 1
   #include "scalarmodels.hh"
   #include "initadvectdiff.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ImplTimeStepper<DgType> ODEType;
#elif PROBLEM == 2
   #include "scalarmodels.hh"
   #include "initburgers.cc"
   typedef U0<GridType> InitialDataType;
   typedef BurgersModel<GridType,InitialDataType > ModelType;
   typedef LLFFlux<ModelType> FluxType;
   //typedef DGLimitedAdvectionOperator<ModelType,LLFFlux,order> DgType;
   typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgType;
   typedef DuneODE::ExplRungeKutta<DgType> ODEType;
#elif PROBLEM == 3
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
#elif PROBLEM == 4
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGLimitedAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DuneODE::ExplTimeStepper<DgType> ODEType;
#endif
