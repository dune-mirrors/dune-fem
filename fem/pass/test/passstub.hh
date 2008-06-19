#ifndef DUNE_PASSSTUB_HH
#define DUNE_PASSSTUB_HH

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>
#include <dune/fem/gridpart/gridpart.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include "dune/fem/pass/pass.hh"
#include "dune/fem/pass/discretemodel.hh"

namespace Dune
{

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

  template< class PreviousPass, int pId = -1 >
  class PassStub
  : public Pass< PassStubTraits< GridType >, PreviousPass, pId >
  {
    typedef PassStub< PreviousPass, pId > ThisType;
    typedef Pass< PassStubTraits< GridType >, PreviousPass, pId > BaseType;

  public:
    typedef PassStubTraits< GridType > PassStubTraitsType;
    typedef PreviousPass PreviousPassType;

    typedef typename BaseType::TotalArgumentType ArgumentType;
    typedef typename PassStubTraitsType::DestinationType DestinationType;

  public:
    explicit PassStub( PreviousPassType &previousPass )
    : BaseType( previousPass )
    {}
      
    virtual void compute(const ArgumentType& arg, DestinationType& dest) const
    {}

    virtual void allocateLocalMemory() {}
  };


  class ProblemStub : public
        DiscreteModelDefaultWithInsideOutSide<PassStubTraits<GridType> > 
  {
  public:  
    typedef PassStubTraits<GridType> Traits;
  };
}

#endif
