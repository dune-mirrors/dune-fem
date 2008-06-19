#ifndef DUNE_DGPASSSTUB_HH
#define DUNE_DGPASSSTUB_HH

#include <string>

#include <dune/grid/io/file/dgfparser/dgfgridtype.hh>

#include "../dgpass.hh"
#include "../discretemodel.hh"
#include "../selection.hh"

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/gridpart/gridpart.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/quadrature/elementquadrature.hh>

namespace Dune {
  
  class DiscreteModelStub;

  struct DGStubTraits {
    enum { dim = GridType :: dimension };
    typedef LeafGridPart<GridType> GridPartType;
    typedef FunctionSpace<double, double, dim , 1> FunctionSpaceType;
    typedef LagrangeDiscreteFunctionSpace<
      FunctionSpaceType, GridPartType, 1> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef ElementQuadrature< GridPartType, 0> VolumeQuadratureType;
    typedef ElementQuadrature< GridPartType, 1>  FaceQuadratureType;
    typedef FieldVector<double, dim > DomainType;
    typedef FieldVector<double, dim > FaceDomainType;
    typedef DiscreteFunctionSpaceType::RangeType RangeType;
    typedef DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    typedef DiscreteModelStub DiscreteModelType;
  };
  
  class DiscreteModelStub : 
    public DiscreteModelDefault<DGStubTraits>
  {
    typedef DiscreteModelDefault<DGStubTraits> BaseType;
  public:
    typedef Selector<0>::Base SelectorType;

    typedef DGStubTraits Traits;

    typedef Traits::FunctionSpaceType FunctionSpaceType;
    typedef Traits::GridPartType GridPartType;
    typedef GridPartType:: GridType GridType;
    typedef Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    // * find common notation
    typedef Traits::SpaceType SpaceType;
    typedef Traits::DestinationType DestinationType;
    typedef Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef Traits::FaceQuadratureType FaceQuadratureType;
    typedef Traits::DiscreteFunctionSpaceType::RangeType RangeType;
    typedef Traits::DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef GridPartType :: IntersectionIteratorType IntersectionIterator;
    typedef GridType::Codim<0>::Entity EntityType;

    typedef BaseType :: MassFactorType MassFactorType; 
 
    typedef FieldVector<double, 3> DomainType;
  
  public:
    bool hasFlux() const { return true; }
    bool hasSource() const { return true; }
    bool hasMass () const { return true; }

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

    template <class ArgumentTuple>
    void mass(const EntityType& en,
              const double time,
              const DomainType& x,
              const ArgumentTuple& u,
              MassFactorType& m)
    { std::cout << "M()" << std::endl; }

  };

} // end namespace Dune

#endif
