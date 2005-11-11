#ifndef DUNE_PASSSTUB_HH
#define DUNE_PASSSTUB_HH

#include <dune/fem/lagrangebase.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/fem/dfadapt.hh>
#include "../pass.hh"

namespace Dune {

  struct PassStubTraits {
    typedef FunctionSpace<double, double, 3, 1> FunctionSpaceType;
    typedef ALU3dGrid<3, 3, tetra> GridType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
    typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
  };

  template <class PreviousPassImp>
  class PassStub : public Pass<PassStubTraits, PreviousPassImp> 
  {
  public:
    typedef Pass<PassStubTraits, PreviousPassImp> BaseType;
    typedef typename BaseType::TotalArgumentType ArgumentType;
    typedef typename PassStubTraits::DestinationType DestinationType;
  public:
    PassStub(PreviousPassImp& prev) :
      BaseType(prev) {}
      
    virtual void compute(const ArgumentType& arg, DestinationType& dest) const
    {}
    virtual void allocateLocalMemory() {}
  private:
  };

}

#endif
