#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include "advectdiff.hh"

#include <dune/fem/solver/rungekutta.hh>
#include <dune/fem/solver/odesolver.hh>
#include <dune/fem/solver/multistep.hh>

using namespace Dune;
// Approximations Ordnung
const int order= POLORDER;
const int rksteps = POLORDER+1; 

// Gitter view Auswahl
//typedef LeafGridPart<GridType> GridPartType;
typedef Dune::HierarchicGridPart<GridType> GridPartType;
//typedef DGAdaptiveLeafGridPart<GridType> GridPartType;

// Modell- und Flussauswahl
// Skalar
#if PROBLEM == 0
   #include "scalarmodels.hh"
   #include "initpulse.cc"
   typedef U0<GridType> InitialDataType;
   //typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   //typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   //typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
   //typedef DGLimitedAdvectionOperator<ModelType,UpwindFlux,order> DgType;
#elif PROBLEM == 1
   #include "scalarmodels.hh"
   #include "initadvectdiff.cc"
   typedef U0<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   typedef LLFFlux<ModelType> FluxType;
   typedef DGAdvectionOperator<ModelType,LLFFlux,order> DgType;
#elif PROBLEM == 2
   #include "scalarmodels.hh"
   #include "initburgers.cc"
   typedef U0<GridType> InitialDataType;
   typedef BurgersModel<GridPartType,InitialDataType > ModelType;
   typedef LLFFlux<ModelType> FluxType;
   //typedef DGLimitedAdvectionOperator<ModelType,LLFFlux,order> DgType;
   //typedef DGAdvectionDiffusionOperator<ModelType,LLFFlux,order> DgType;
   typedef DGAdvectionOperator<ModelType,LLFFlux,order> DgType;
#elif PROBLEM == 3
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionDiffusionOperator<ModelType,UpwindFlux,order> DgType;
#elif PROBLEM == 4
#include "scalarmodels.hh"
#include "initadvectdiff.cc"
   typedef U0Disc<GridType> InitialDataType;
   typedef AdvectionDiffusionModel<GridPartType,InitialDataType> ModelType;
   // typedef LLFFlux<ModelType> FluxType;
   typedef UpwindFlux<ModelType> FluxType;
   typedef DGAdvectionOperator<ModelType,UpwindFlux,order> DgType;
#elif PROBLEM == 5
#include "euler_mhd/eulermodel.hh"
   typedef U0RotatingCone InitialDataType;
   typedef EulerModel<GridPartType,InitialDataType> ModelType;
   typedef HLLNumFlux<ModelType> FluxType;
   // typedef DWNumFlux<ModelType> FluxType;
   typedef DGAdvectionOperator<ModelType,HLLNumFlux,order> DgType;
#elif PROBLEM == 6
#include "euler_mhd/eulermodel.hh"
   typedef U0RP InitialDataType;
   typedef EulerModel<GridPartType,InitialDataType> ModelType;
   typedef HLLNumFlux<ModelType> FluxType;
   // typedef DWNumFlux<ModelType> FluxType;
   typedef DGAdvectionOperator<ModelType,HLLNumFlux,order> DgType;
#endif
   
typedef DuneODE::ExplicitOdeSolver<DgType::DestinationType> ODEType;
// typedef DuneODE::ImplicitOdeSolver<DgType::DestinationType> ODEType;

typedef Tuple< DgType :: DestinationType * > IOTupleType;
