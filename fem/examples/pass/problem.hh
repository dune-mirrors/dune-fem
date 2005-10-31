#ifndef DUNE_EXAMPLEPROBLEM_HH
#define DUNE_EXAMPLEPROBLEM_HH

#include "../../pass/problem.hh"
#include "../../pass/selection.hh"

// Dune includes
#include <dune/common/utility.hh>
#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/discretefunction/adaptivefunction.hh>
#include <dune/grid/alu3dgrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

using namespace Dune;

namespace Dune {

  struct TransportDiffusionTraits1 
  {
    enum { dimRange = 3 };
    enum { dimDomain = 3 };
    enum { polOrd = 0 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef FunctionSpace<
      double, double, dimDomain, 1> SingleFunctionSpaceType;
    typedef ALU3dGrid<dimDomain, dimDomain, tetra> GridType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef DiscontinuousGalerkinSpace<
      SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    // * Need to do: adapt quadrature order
    typedef FixedOrderQuad<
      double, FieldVector<double, 3>, 5> VolumeQuadratureType;
    typedef FixedOrderQuad<
      double, FieldVector<double, 2>, 5> FaceQuadratureType;
  };

  class TransportDiffusionProblem1 :
    public ProblemDefault<
    TransportDiffusionProblem1,
    TransportDiffusionTraits1::FunctionSpaceType> {
  public:
    typedef TransportDiffusionTraits1 Traits;

    typedef Selector<0>::Base SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

  public:
    TransportDiffusionProblem1(double epsilon) :
      epsilon_(epsilon)
    {}

    bool hasSource() const { return true; }
    bool hasFlux() const { return false; }
    
    template <
      class Entity, class ArgumentTuple, class JacobianTuple, class ResultType>
    void source(Entity& en, double time, const DomainType& x, 
                const ArgumentTuple& u, const JacobianTuple& gradU, 
                ResultType& s) {
      typedef typename ElementType<0, JacobianTuple>::Type GradientType;
      s = (Element<0>::get(gradU))[0];
      s *= epsilon_;
    }

  private:
    double epsilon_;
  };

  struct TransportDiffusionTraits2 
  {
    enum { dimRange = 1 };
    enum { dimDomain = 3 };
    enum { polOrd = 0 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef ALU3dGrid<dimDomain, dimDomain, tetra> GridType;
    typedef LeafGridPart<GridType> GridPartType;
    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    // * Need to do: adapt quadrature order
    typedef FixedOrderQuad<
      double, FieldVector<double, 3>, 5> VolumeQuadratureType;
    typedef FixedOrderQuad<
      double, FieldVector<double, 2>, 5> FaceQuadratureType;

  };

  class TransportDiffusionProblem2 : 
    public ProblemDefault<
    TransportDiffusionProblem2,
    TransportDiffusionTraits2::FunctionSpaceType> {
  public:
    typedef TransportDiffusionTraits2 Traits;
    typedef Selector<0, 1>::Base SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
  public:
    TransportDiffusionProblem2(const DomainType& velocity) :
      velocity_(velocity)
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <
      class IntersectionIterator, class ArgumentTuple, class ResultType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         ResultType& gLeft,
                         ResultType& gRight)
    { 
      double upwind = it.unitOuterNormal(x)*velocity_;
      if (upwind > 0) {
        analyticalFlux(*it.inside(), time, it.intersectionGlobal().global(x),
                       uLeft, gLeft);
        gRight = -gLeft;
      } else {
        analyticalFlux(*it.outside(), time, it.intersectionGlobal().global(x),
                       uRight, gLeft);
        gRight = -gLeft;
      }
      return upwind;
    }

    template <class Entity, class ArgumentTuple, class ResultType>
    void analyticalFlux(Entity& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, ResultType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type W1Type;

      const UType& argU = Element<0>::get(u);
      const W1Type& argW1 = Element<1>::get(u);
      //f = argU*velocity_ - argW1;
      //f[0] = velocity_;
      //f *= argU;
      //f[0] -= argW1;
    }

  private:
    DomainType velocity_;

  };

}

#endif
