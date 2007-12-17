#ifndef DUNE_PASSSTUB_HH
#define DUNE_PASSSTUB_HH

#include <dune/fem/space/lagrangespace.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include "../pass.hh"
#include "../discretemodel.hh"

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

namespace Dune {

  class ProblemStub;

  template <class GridImp> 
  struct PassStubTraits {
    typedef GridImp GridType;
    typedef PassStubTraits<GridType> Traits;
    typedef ProblemStub DiscreteModelType;
    typedef FunctionSpace<double, double, GridType :: dimension, 1> FunctionSpaceType;
    typedef typename FunctionSpaceType :: RangeType  RangeType;
    typedef typename FunctionSpaceType :: DomainType DomainType;
    typedef typename FunctionSpaceType :: JacobianRangeType JacobianRangeType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef ElementQuadrature<GridPartType,0> VolumeQuadratureType;
    typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
  };

  template <class PreviousPassImp>
  class PassStub : public Pass<PassStubTraits<GridType> , PreviousPassImp> 
  {
  public:
    typedef PassStubTraits<GridType> PassStubTraitsType;
    typedef Pass<PassStubTraitsType , PreviousPassImp> BaseType;
    typedef typename BaseType::TotalArgumentType ArgumentType;
    typedef typename PassStubTraitsType::DestinationType DestinationType;
  public:
    PassStub(PreviousPassImp& prev) :
      BaseType(prev) {}
      
    virtual void compute(const ArgumentType& arg, DestinationType& dest) const
    {}
    virtual void allocateLocalMemory() {}
  private:
  };


  class ProblemStub : public
        DiscreteModelDefaultWithInsideOutSide<PassStubTraits<GridType> > 
  {
  public:  
    typedef PassStubTraits<GridType> Traits;
  };
}

#endif
