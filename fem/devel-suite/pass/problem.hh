#ifndef DUNE_EXAMPLEPROBLEM_HH
#define DUNE_EXAMPLEPROBLEM_HH

#include "../../pass/discretemodel.hh"
#include "../../pass/selection.hh"

// Dune includes
#include <dune/common/utility.hh>
#include "../../space/dgspace.hh"
#include <dune/fem/space/combinedspace.hh>
#include <dune/fem/dfadapt.hh>
#include <dune/fem/function/adaptivefunction.hh>

#if HAVE_ALUGRID
#include <dune/grid/alu3dgrid.hh>
#endif

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#include <dune/grid/sgrid.hh>
#include <dune/fem/gridpart/gridpart.hh>
//#include <dune/quadrature/fixedorder.hh>
#include "../../quadrature/quadrature.hh"
#include "../../quadrature/cachequad.hh"

namespace Dune {

  class TransportDiffusionDiscreteModel1;
  class TransportDiffusionDiscreteModel2;

  struct TransportDiffusionTraits1 
  {
    enum { dimRange = 2 };
    enum { dimDomain = 2 };
    enum { polOrd = 2 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;
    typedef FunctionSpace<
      double, double, dimDomain, 1> SingleFunctionSpaceType;

    //typedef AlbertaGrid<dimDomain, dimDomain> GridType;
    typedef SGrid<dimDomain, dimDomain> GridType;
    typedef LeafGridPart<GridType> GridPartType;
    //typedef DiscontinuousGalerkinSpace<
    //  SingleFunctionSpaceType, GridPartType, polOrd> SingleSpaceType;
    //typedef CombinedSpace<SingleSpaceType, dimRange> DiscreteFunctionSpaceType;
    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef FunctionSpaceType::DomainType DomainType;
    typedef FunctionSpaceType::RangeType RangeType;
    typedef FunctionSpaceType::JacobianRangeType JacobianRangeType;

    // * Need to do: adapt quadrature order
    /*
    typedef Quadrature<double, dimDomain> VolumeQuadratureType;
    typedef Quadrature<double, dimDomain-1> FaceQuadratureType;
    */
    /*
    typedef CacheQuadrature<GridType, 0> VolumeQuadratureType;
    typedef CacheQuadrature<GridType, 1> FaceQuadratureType;
    */
    typedef CachingQuadrature<GridType, 0> VolumeQuadratureType;
    typedef CachingQuadrature<GridType, 1> FaceQuadratureType;    

    typedef TransportDiffusionDiscreteModel1 DiscreteModelType;
  };

  class TransportDiffusionDiscreteModel1 :
    public DiscreteModelDefault<TransportDiffusionTraits1> {
  public:
    typedef TransportDiffusionTraits1 Traits;

    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef Traits::RangeType RangeType;
    typedef Traits::GridType GridType;
    typedef GridType::Codim<0>::Entity EntityType;

  public:
    TransportDiffusionDiscreteModel1() :
      identity_(0.0)
    {
      identity_[0][0] = identity_[1][1] = 1.0;
    }

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }
    
    template <class ArgumentTuple, class JacobianTuple>
    void source(EntityType& en, double time, const DomainType& x, 
                const ArgumentTuple& u, const JacobianTuple& gradU, 
                RangeType& s) {
      typedef typename ElementType<0, JacobianTuple>::Type GradientType;
      //std::cout << "x == " << en.geometry().global(x) << ": " 
      //          << Element<0>::get(u)[0] << ", "
      //          << Element<0>::get(gradU)[0] << std::endl;
    }

    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight)
    { 
      // central differences
      //std::cout << "x == " << x << ": " << Element<0>::get(uLeft)[0]
      //           << " ?= " << Element<0>::get(uRight)[0] << std::endl;

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
    enum { polOrd = 2 };

    typedef FunctionSpace<
      double, double, dimDomain, dimRange> FunctionSpaceType;

    //typedef AlbertaGrid<dimDomain, dimDomain> GridType;
    typedef SGrid<dimDomain, dimDomain> GridType;
    //typedef LeafGridPart<GridType> GridPartType;
    typedef DefaultGridIndexSet<GridType, LevelIndex> IndexSetType;
    typedef DefaultGridPart<GridType, IndexSetType> GridPartType;
    typedef DiscontinuousGalerkinSpace<
      FunctionSpaceType, GridPartType, polOrd> DiscreteFunctionSpaceType;
    typedef AdaptiveDiscreteFunction<DiscreteFunctionSpaceType> DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef FunctionSpaceType::DomainType DomainType;
    typedef FunctionSpaceType::RangeType RangeType;
    typedef FunctionSpaceType::JacobianRangeType JacobianRangeType;

    // * Need to do: adapt quadrature order
    /*
    typedef Quadrature<double, dimDomain> VolumeQuadratureType;
    typedef Quadrature<double, dimDomain-1> FaceQuadratureType;
    */
    /*
    typedef CacheQuadrature<GridType, 0> VolumeQuadratureType;
    typedef CacheQuadrature<GridType, 1> FaceQuadratureType;
    */
    typedef CachingQuadrature<GridType, 0> VolumeQuadratureType;
    typedef CachingQuadrature<GridType, 1> FaceQuadratureType;

    typedef TransportDiffusionDiscreteModel2 DiscreteModelType;
  };

  class TransportDiffusionDiscreteModel2 : 
    public DiscreteModelDefault<TransportDiffusionTraits2> {
  public:
    typedef TransportDiffusionTraits2 Traits;
    typedef Selector<0, 1> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;

    typedef Traits::GridType GridType;
    typedef Traits::DiscreteFunctionType::JacobianRangeType JacobianRangeType;
    typedef GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef GridType::Codim<0>::Entity EntityType;

  public:
    TransportDiffusionDiscreteModel2(const DomainType& velocity, double epsilon) :
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
      gLeft *= upwind;
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
      f *= argU[0];
      f[0] += argW1;
 
    }

  private:
    DomainType velocity_;
    double epsilon_;
  };

}

#endif
