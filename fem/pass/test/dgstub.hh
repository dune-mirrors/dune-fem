#ifndef DUNE_DGPASSSTUB_HH
#define DUNE_DGPASSSTUB_HH

#include <string>

#include "../dgpass.hh"
#include "../discretemodel.hh"
#include "../selection.hh"

#include <dune/fem/lagrangebase.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/quadrature/fixedorder.hh>

namespace Dune {
  
  class DiscreteModelStub;

  struct DGStubTraits {
    typedef FunctionSpace<double, double, 3, 1> FunctionSpaceType;
    typedef ALU3dGrid<3, 3, tetra> GridType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    //typedef DFAdapt<DiscreteFunctionSpaceType> DestinationType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef FixedOrderQuad<double, FieldVector<double, 3>, 1> VolumeQuadratureType;
    typedef FixedOrderQuad<double, FieldVector<double, 2>, 1> FaceQuadratureType;
    typedef FieldVector<double, 3> DomainType;
    typedef FieldVector<double, 2> FaceDomainType;
    typedef DiscreteFunctionSpaceType::RangeType RangeType;
    typedef DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscreteModelStub DiscreteModelType;
  };
  
  class DiscreteModelStub : 
    public DiscreteModelDefault<DGStubTraits>
  {
  public:
    typedef Selector<0>::Base SelectorType;

    typedef DGStubTraits Traits;

    typedef Traits::FunctionSpaceType FunctionSpaceType;
    typedef Traits::GridType GridType;
    typedef Traits::GridPartType GridPartType;
    typedef Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    // * find common notation
    typedef Traits::SpaceType SpaceType;
    typedef Traits::DestinationType DestinationType;
    typedef Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef Traits::FaceQuadratureType FaceQuadratureType;
    typedef Traits::DiscreteFunctionSpaceType::RangeType RangeType;
    typedef Traits::DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef Traits::GridType GridType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef GridType::Codim<0>::Entity EntityType;
 
    typedef FieldVector<double, 3> DomainType;
  
  public:
    bool hasFlux() const { return true; }
    bool hasSource() const { return true; }

    template <class ArgumentTuple, class FaceDomainType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      std::cout << "numericalFlux()" << std::endl; 
      return 0.0;
    }

    template <class ArgumentTuple, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& boundaryFlux)
    {
      std::cout << "boundaryFlux()" << std::endl;
      return 0.0;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f) 
    { std::cout << "analyticalFlux()" << std::endl; }

    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, 
                double time, const DomainType& x,
                const ArgumentTuple& u, 
                const JacobianTuple& jac, 
                RangeType& s)
    { std::cout << "S()" << std::endl; }
  };

} // end namespace Dune

#endif
