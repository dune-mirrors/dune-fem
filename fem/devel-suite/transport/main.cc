// Global defines

#define SGRID 0
#define AGRID 1

#include "../../config.h"

#if !HAVE_ALBERTA 
#undef SGRID 
#undef AGRID 
#define SGRID 1
#define AGRID 0
#endif 

#if SGRID
#include <dune/grid/sgrid.hh>
#define GridName SGrid 
#endif 

#if AGRID
#include <dune/grid/albertagrid.hh>
#define GridName AlbertaGrid 
#endif

#define DIM 2
#define DIM_OF_WORLD 2


// Include system headers
#include <iostream>
#include <strstream>
#include <string>

// Include Dune headers
#include <dune/fem/discreteoperatorimp.hh>
#include <dune/fem/lagrangebase.hh>
#include <dune/grid/common/leafindexset.hh>
#include <dune/fem/inverseoperators.hh>
#include <dune/fem/discfuncarray.hh>

#include <dune/io/file/grapedataio.hh>

// Include local headers
#include "description.hh"
#include "../../misc/advection.hh"
#include "riemann.cc"
#include "../../misc/inverseoperatorfactory.hh"
#include "../../misc/timestepping.hh"

#include <dune/common/stdstreams.cc>

using namespace Dune;
using namespace Adi;

const int dimRange = 3;

// Typedefs
typedef Descr<double,
              DIM, 
              DIM_OF_WORLD, 
              GridName, 
              1, 
              RiemannProblem, 
              0> DescrType;
typedef DescrType::DiscreteFunctionType DiscreteFunction;
typedef TransportDescription<DescrType> ProblemType;
typedef DescrType::GridType GridType;

void printData(double time,
               int timestep,
               ProblemType& pd,
               ExplicitEuler<ProblemType>& loop) {
  std::stringstream stream0;
  stream0 << "ascii" << timestep;

  GrapeDataIO<GridType> dataio;
  dataio.writeGrid(pd.grid(), xdr, "grid", time, timestep);
  
   //  sol0.write_xdr(stream2.str().c_str());
  //sol1.write_xdr(stream3.str().c_str());
  dataio.writeData(loop.solution(), xdr, "vec", timestep);
  loop.solution().write_ascii(stream0.str().c_str());
 }

int main() {
  // Parameter file: parameter
  std::string param("parameter");
  double dt = 0.01;

  // Generate problem object from parameter file
  ProblemType problem(param.c_str());

  // Create timestepping
  CGInverseOperatorFactory<DiscreteFunction> factory(1E-6, 1E-10, 100000, 0);
  ExplicitEuler<ProblemType> loop(problem, factory, dt);

  printData(loop.time(), 0, problem, loop);
  loop.solve();
  printData(loop.time(), 1, problem, loop);
  loop.solveUntil(problem.endTime());
  printData(loop.time(), 2, problem, loop);

}
