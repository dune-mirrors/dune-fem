#ifndef DUNE_EXAMPLEPROBLEM_HH
#define DUNE_EXAMPLEPROBLEM_HH

#include "../../pass/problem.hh"
#include <dune/common/utility.hh>

using namespace Dune;

namespace Adi {

  class TransportDiffusionProblem1 
    public ProblemDefault {
  public:
    TransportDiffusionProblem1(double epsilon) :
      epsilon_(epsilon)
    {}

    bool hasSource() const { return true; }
    bool hasFlux() const { return false; }
    
    template <
      class Entity, class ArgumentTuple, class JacobianTuple, class ResultType
    void source(Entity& en, double time, const DomainType& x, 
                const ArgumentTuple& u, const JacobianTuple& gradU, 
                ResultType& s) {
      typedef typename ElementType<0, JacobianTuple>::Type GradientType;
      s = epsilon_*Element<0>::get(gradU);
    }

  private:
    double epsilon_;
  };

  class TransportDiffusionProblem2 : 
    public ProblemDefault {
  public:
    TransportDiffusionProblem2(const DomainType& velocity) :
      velocity_(velocity)
    {}

    bool hasSource() const { return false; }
    bool hasFlux() const { return true; }

    template <
      class IntersectionIterator, class ArgumentTuple, class ResultType>
    double numericalFlux(IntersectionIterator& it,
                         double time, const DomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         ResultType& gLeft,
                         ResultType& gRight)
    { 
      double upwind = it.outerNormal()*velocity_;
      if (upwind > 0) {
        analyticalFlux(*it.inside(), time, x, uLeft, gLeft);
        gRight = -gLeft;
      } else {
        analyticalFlux(*it.outside(), time, x, uRight, gLeft);
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

      f = velocity_*Element<0>::get(u) - Element<1>::get(u);
    }

  private:
    DomainType velocity_;

  };

}

#endif
