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
#include <dune/grid/albertagrid.hh>
#include <dune/grid/common/gridpart.hh>
#include <dune/quadrature/fixedorder.hh>

using namespace Dune;

namespace Dune {

  class TransportDiffusionProblem1;
  class TransportDiffusionProblem2;

  struct TransportDiffusionTraits1 
  {
    enum { dimRange = 2 };
    enum { dimDomain = 2 };
    enum { polOrd = 3 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef FunctionSpace<
      double, double, dimDomain, 1> SingleFunctionSpaceType;

    //    typedef AlbertaGrid<dimDomain, dimDomain> GridType;
    typedef SGrid<dimDomain, dimDomain> GridType;
    //typedef LeafGridPart<GridType> GridPartType;
    typedef DefaultGridIndexSet<GridType, LevelIndex> IndexSetType;
    typedef DefaultGridPart<GridType, IndexSetType> GridPartType;
    typedef DiscontinuousGalerkinSpace<
      SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef FunctionSpaceType::DomainType DomainType;
    typedef FunctionSpaceType::RangeType RangeType;
    typedef FunctionSpaceType::JacobianRangeType JacobianRangeType;

    // * Need to do: adapt quadrature order
    typedef FixedOrderQuad<
      double, FieldVector<double, dimDomain>, 5> VolumeQuadratureType;
    typedef FixedOrderQuad<
      double, FieldVector<double, dimDomain-1>, 5> FaceQuadratureType;

    typedef TransportDiffusionProblem1 ProblemType;
  };

  class TransportDiffusionProblem1 :
    public ProblemDefault<TransportDiffusionTraits1> {
  public:
    typedef TransportDiffusionTraits1 Traits;

    typedef Selector<0>::Base SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef Traits::RangeType RangeType;
    typedef Traits::GridType GridType;
    typedef GridType::Codim<0>::Entity EntityType;

  public:
    TransportDiffusionProblem1() :
      identity_(0.0)
    {
      identity_[0][0] = identity_[1][1] = 1.0;
    }

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }
    
    /*
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, double time, const DomainType& x, 
                const ArgumentTuple& u, const JacobianTuple& gradU, 
                RangeType& s) {
      typedef typename ElementType<0, JacobianTuple>::Type GradientType;
      s = (Element<0>::get(gradU))[0];
    }
    */

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      // central differences
      const DomainType normal = it.integrationOuterNormal(x);
      double uMean = 
        0.5*( (Element<0>::get(uLeft))[0] + (Element<0>::get(uRight))[0]);
   
      gLeft = normal;
      gRight = normal;
      gLeft *= uMean;
      gRight *= uMean;
      
      return gLeft*gLeft;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      JacobianRangeType anaFlux(0.0);
      DomainType normal = it.integrationOuterNormal(x);
      analyticalFlux(*it.inside(), time, it.intersectionSelfLocal().global(x),
                     uLeft, anaFlux);
      gLeft *= 0.0;
      anaFlux.umv(normal, gLeft);
      return 1.0; // dummy, need to change that
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      f = identity_;
      f *= (Element<0>::get(u))[0];
    }

  private:
    JacobianRangeType identity_;
  };

  struct TransportDiffusionTraits2 
  {
    enum { dimRange = 1 };
    enum { dimDomain = 2 };
    enum { polOrd = 3 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;

    //    typedef AlbertaGrid<dimDomain, dimDomain> GridType;
    typedef SGrid<dimDomain, dimDomain> GridType;
    //typedef LeafGridPart<GridType> GridPartType;
    typedef DefaultGridIndexSet<GridType, LevelIndex> IndexSetType;
    typedef DefaultGridPart<GridType, IndexSetType> GridPartType;
    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    typedef DiscreteFunctionSpaceType SpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef FunctionSpaceType::DomainType DomainType;
    typedef FunctionSpaceType::RangeType RangeType;
    typedef FunctionSpaceType::JacobianRangeType JacobianRangeType;

    // * Need to do: adapt quadrature order
    typedef FixedOrderQuad<
      double, FieldVector<double, dimDomain>, 5> VolumeQuadratureType;
    typedef FixedOrderQuad<
      double, FieldVector<double, dimDomain-1>, 5> FaceQuadratureType;

    typedef TransportDiffusionProblem2 ProblemType;
  };

  class TransportDiffusionProblem2 : 
    public ProblemDefault<TransportDiffusionTraits2> {
  public:
    typedef TransportDiffusionTraits2 Traits;
    typedef Selector<0, 1>::Base SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef Traits::GridType GridType;
    typedef Traits::DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef GridType::Codim<0>::Entity EntityType;

  public:
    TransportDiffusionProblem2(const DomainType& velocity, double epsilon) :
      velocity_(velocity),
      epsilon_(epsilon)
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      DomainType normal = it.integrationOuterNormal(x);
      double upwind = normal*velocity_;
      
      // transport contribution
      if (upwind > 0) {
        gLeft = Element<0>::get(uLeft);
      }
      else {
        gLeft = Element<0>::get(uRight);
      }
      gLeft *= -upwind;
      gRight = gLeft;

      // diffusion contribution
      DomainType vMean = 
        Element<1>::get(uLeft) + Element<1>::get(uRight);
      vMean *= 0.5*epsilon_;

      gLeft[0] += normal*vMean;
      gRight[0] += normal*vMean;
    
      /* old version
      const DomainType normal = it.integrationOuterNormal(x);
      double upwind = normal*velocity_;
      JacobianRangeType anaFlux(0.0);
      
      if (upwind > 0) {
        analyticalFlux(*it.inside(), time, 
                       it.intersectionSelfLocal().global(x),
                       uLeft, anaFlux);
      } else {
        analyticalFlux(*it.outside(), time, 
                       it.intersectionSelfLocal().global(x),
                       uRight, anaFlux);
      }
      gLeft *= 0.0;
      anaFlux.umv(normal, gLeft);
      gRight = gLeft;
      /*
      std::cout << "numericalFlux:\n";
      std::cout << Element<0>::get(uLeft) << ", " << Element<0>::get(uRight) << std::endl;
      std::cout << anaFlux << std::endl;
      std::cout << normal << std::endl;
      std::cout << "result = " << gLeft << std::endl;
      */
      return upwind;
    }

    template <class ArgumentTuple>
    double boundaryFlux(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft,
                        RangeType& gLeft)
    {
      // Just make an outflow boundary for now
      JacobianRangeType anaFlux;
      DomainType normal = it.integrationOuterNormal(x);
      analyticalFlux(*it.inside(), time, it.intersectionSelfLocal().global(x),
                     uLeft, anaFlux);
      
      gLeft *= 0.0;
      anaFlux.umv(normal, gLeft);

      return normal*velocity_;
    }

    template <class ArgumentTuple>
    void analyticalFlux(EntityType& en,
                        double time, const DomainType& x,
                        const ArgumentTuple& u, JacobianRangeType& f)
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      typedef typename ElementType<1, ArgumentTuple>::Type W1Type;

      const UType& argU = Element<0>::get(u);
      W1Type argW1 = Element<1>::get(u);
      argW1 *= epsilon_;

      //f = argU*velocity_ + argW1;
      f[0] = velocity_;
      f *= -argU[0];
      f[0] += argW1;
 
    }

  private:
    DomainType velocity_;
    double epsilon_;
  };

}

#endif
