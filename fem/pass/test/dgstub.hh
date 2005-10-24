#ifndef DUNE_DGPASSSTUB_HH
#define DUNE_DGPASSSTUB_HH

#include <string>

#include "../dgpass.hh"
#include "../problem.hh"
#include "../selection.hh"

#include <dune/fem/lagrangebase.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/quadrature/fixedorder.hh>

namespace Dune {

 struct DGStubTraits {
   typedef FunctionSpace<double, double, 3, 1> FunctionSpaceType;
   typedef ALU3dGrid<3, 3, tetra> GridType;
   typedef LeafGridPart<GridType> GridPartType;
   typedef LagrangeDiscreteFunctionSpace<
     FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
   typedef DiscreteFunctionSpaceType SpaceType;
   typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
   typedef FixedOrderQuad<double, FieldVector<double, 3>, 1> VolumeQuadratureType;
   typedef FixedOrderQuad<double, FieldVector<double, 2>, 1> FaceQuadratureType;
  };
  
  class ProblemStub : 
    public ProblemDefault<ProblemStub, DGStubTraits::FunctionSpaceType>
  {
  public:
    typedef Selector<0>::Base SelectorType;

  public:
    bool hasFlux() const { return true; }
    bool hasNonConservativeTerm() const { return true; }
    bool hasSource() const { return true; }

    void f() { std::cout << "f()" << std::endl; }
    void g() { std::cout << "g()" << std::endl; }
    void S() { std::cout << "S()" << std::endl; }
  };

} // end namespace Dune

#endif
