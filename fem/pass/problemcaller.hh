#ifndef DUNE_PROBLEMCALLER_HH
#define DUNE_PROBLEMCALLER_HH

#include <utility>

#include "caller.hh"
#include "problem.hh"

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class ProblemImp, class ArgumentImp, class SelectorImp>
  class ProblemCaller {
  public:
    typedef ProblemImp ProblemType;
    typedef ArgumentImp TotalArgumentType;
    typedef SelectorImp SelectorType;

    typedef Filter<TotalArgumentType, SelectorType> FilterType;
    typedef typename FilterType::ResultType DiscreteFunctionTupleType;
    typedef typename LocalFunctionCreator<
      DiscreteFunctionTupleType>::ResultType LocalFunctionTupleType;
    typedef typename RangeVectorCreator<
      LocalFunctionTupleType>::ResultType RangeTupleType;
    //typedef typename Caller<RangeTupleType> CallerType;

    typedef LocalFunctionCreator<DiscreteFunctionTupleType> LFCreator;
    typedef RangeVectorCreator<LocalFunctionTupleType> RangeCreator;

    typedef std::pair<LocalFunctionTupleType, RangeTupleType> ValuePair;
  public:
    ProblemCaller(ProblemType& problem, TotalArgumentType& arg) :
      problem_(problem),
      arg_(&arg),
      discreteFunctions_(FilterType::apply(arg)),
      valuesEn_(LFCreator::apply(discreteFunctions_), RangeCreator::apply()),
      valuesNeigh_(LFCreator::apply(discreteFunctions_), RangeCreator::apply())
    {}

    void setArgument(TotalArgumentType& arg) {      
      if (*arg_ != arg) {

        // Set pointer
        arg_ = &arg;
      
        // Filter the argument
        discreteFunctions_ = FilterType::apply(arg);

        // Build up new local function tuples
        valuesEn_.first = LFCreator::apply(discreteFunctions_);
        valuesNeigh_.first = LFCreator::apply(discreteFunctions_);
      }
    }

    template <class Entity>
    void setEntity(Entity& en) {
      setter(en, valuesEn_.first);
    }

    template <class Entity>
    void setNeighbor(Entity& en) {
      setter(en, valuesNeigh_.first);
    }

    // Here, the interface of the problem is replicated and the Caller
    // is used to do the actual work
    template <class Entity, class DomainType, class ResultType>
    void analyticalFlux(Entity& en, const DomainType& x, ResultType& res) {
      problem_.analyticalFlux();
      //CallerType::analyticalFlux(problem_, en, x, valuesEn_.second, res);
    }
    
    template <class Entity, class QuadratureType, class ResultType>
    void analyticalFlux(Entity& en, QuadratureType& quad, int quadPoint,
                        ResultType& res) {
      problem_.analyticalFlux();
      //CallerType::analyticalFlux(problem_, en, quad, quadPoint, valuesEn_.second, res);
    }

    template <class IntersectionIterator, class DomainType, class ResultType>
    void numericalFlux(IntersectionIterator& nit,
                       const DomainType& x,
                       ResultType& res_en, ResultType& res_neigh) {
      problem_.numericalFlux();
      //CallerType::numericalFlux(problem_, nit, x, 
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }

    template <class IntersectionIterator, class QuadratureType, class ResultType>
    void numericalFlux(IntersectionIterator& nit,
                       QuadratureType& quad,
                       int quadPoint,
                       ResultType& res_en, ResultType& res_neigh) {
      problem_.numericalFlux();
      //CallerType::numericalFlux(problem_, nit, quad, quadPoint,
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }

    template <class Entity, class DomainType, class ResultType>
    void source(Entity& en, const DomainType& x, ResultType& res) {
      problem_.source();
      //CallerType::source(problem_, en, x, valuesEn_.second, res);
    }

    template <class Entity, class QuadratureType, class ResultType>
    void source(Entity& en, QuadratureType& quad, int quadPoint, ResultType& res) {
      problem_.source();
      //CallerType::source(problem_, en, quad, quadPoint, valuesEn_.second, res);
    }

  private:
    template <class Entity>
    void setter(Entity& en, LocalFunctionTupleType& tuple) {
      ForEachValue<LocalFunctionTupleType> forEach(tuple);
      forEach.apply(LocalFunctionSetter<Entity>(en));
    }

    template <class Entity, class DomainType>
    void evaluateLocal(Entity& en, const DomainType& x, ValuePair& p) {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(p.first,
                                                              p.second);

      forEach.apply(LocalFunctionEvaluateLocal<Entity, DomainType>(en, x));
    }

    template <class Entity, class QuadratureType>
    void evaluateQuad(Entity& en, QuadratureType& quad, int quadPoint, 
                      ValuePair& p) {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(p.first,
                                                              p.second);

      forEach.apply(LocalFunctionEvaluateQuad<
                    Entity, QuadratureType>(en, quad, quadPoint));
    }

  private:
    ProblemCaller(const ProblemCaller&);
    ProblemCaller& operator=(const ProblemCaller&);

  private:
    ProblemType& problem_;
    TotalArgumentType* arg_;

    DiscreteFunctionTupleType discreteFunctions_;
    //LocalFunctionTupleType localFunctionsEn_;
    //LocalFunctionTupleType localFunctionsNeigh_;
    //RangeTupleType rangesEn_;
    //RangeTupleType rangesNeigh_;
    ValuePair valuesEn_;
    ValuePair valuesNeigh_;
  };

}

#endif
