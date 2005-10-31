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
    typedef typename Creator<
      RangeTypeEvaluator,
      LocalFunctionTupleType>::ResultType RangeTupleType;
    typedef typename Creator<
      JacobianRangeTypeEvaluator,
      LocalFunctionTupleType>::ResultType JacobianRangeTupleType;
    //typedef typename Caller<RangeTupleType> CallerType;

    typedef LocalFunctionCreator<DiscreteFunctionTupleType> LFCreator;
    typedef Creator<
      RangeTypeEvaluator, LocalFunctionTupleType> RangeCreator;
    typedef Creator<
      JacobianRangeTypeEvaluator, LocalFunctionTupleType> JacobianCreator;

    typedef std::pair<LocalFunctionTupleType, RangeTupleType> ValuePair;
  public:
    ProblemCaller(ProblemType& problem, TotalArgumentType& arg) :
      problem_(problem),
      arg_(&arg),
      discreteFunctions_(FilterType::apply(arg)),
      valuesEn_(LFCreator::apply(discreteFunctions_), RangeCreator::apply()),
      valuesNeigh_(LFCreator::apply(discreteFunctions_),RangeCreator::apply()),
      jacobians_(JacobianCreator::apply())
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
      evaluateLocal(en, x, valuesEn_);
      problem_.analyticalFlux(en, 0.0, x, valuesEn_.second, res);
      //CallerType::analyticalFlux(problem_, en, x, valuesEn_.second, res);
    }
    
    template <class Entity, class QuadratureType, class ResultType>
    void analyticalFlux(Entity& en, QuadratureType& quad, int quadPoint,
                        ResultType& res) {
      // * temporary
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(valuesEn_.first,
                                                        valuesEn_.second);
      LocalFunctionEvaluateQuad<Entity, QuadratureType> eval(en, quad, quadPoint);
      forEach.apply(eval);


      //evaluateQuad(en, quad, quadPoint, valuesEn_);
      problem_.analyticalFlux(en, 0.0, quad.point(quadPoint), valuesEn_.second, res);
      //CallerType::analyticalFlux(problem_, en, quad, quadPoint, valuesEn_.second, res);
    }

    template <class IntersectionIterator, class DomainType, class ResultType>
    double numericalFlux(IntersectionIterator& nit,
                         const DomainType& x,
                         ResultType& resEn, ResultType& resNeigh) {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      const Geometry& neighLocal = nit.intersectionNeighborLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(x),
                    valuesEn_);
      evaluateLocal(*nit.outside(), neighLocal.global(x),
                    valuesNeigh_);
      return problem_.numericalFlux(nit, 0.0, x,
                                    valuesEn_.second, valuesNeigh_.second,
                                    resEn, resNeigh);
      //CallerType::numericalFlux(problem_, nit, x, 
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }

    // Ensure: entities set correctly before call
    template <
      class IntersectionIterator, class QuadratureType, class ResultType>
    double numericalFlux(IntersectionIterator& nit,
                         QuadratureType& quad,
                         int quadPoint,
                         ResultType& resEn, ResultType& resNeigh) {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      const Geometry& neighLocal = nit.intersectionNeighborLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(quad.point(quadPoint)),
                    valuesEn_);
      evaluateLocal(*nit.outside(), neighLocal.global(quad.point(quadPoint)),
                    valuesNeigh_);
      return problem_.numericalFlux(nit, 0.0, 
                                    quad.point(quadPoint),
                                    valuesEn_.second, valuesNeigh_.second,
                                    resEn, resNeigh);
      //CallerType::numericalFlux(problem_, nit, quad, quadPoint,
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }

    // * Problem: bval should be a tuple as well, but how to you get it?
    /*
    template <class IntersectionIterator, class DomainType, class RangeType, class ResultType>
    double boundaryFlux(IntersectionIterator& nit,
                        const DomainType& x,
                        const RangeType& bval,
                        ResultType& resEn) {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      ResultType dummy(0.0);
      const Geometry& global = nit.intersectionGlobal();
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(x),
                    valuesEn_);
      return problem_.numericalFlux(nit, 0.0,
                                    global.global(x),
                                    valuesEn_.second, bval,
                                    resEn, dummy);
    }

    template <
      class IntersectionIterator, class QuadratureType, class RangeType, class ResultType>
    double boundaryFlux(IntersectionIterator& nit,
                        QuadratureType& quad,
                        int quadPoint,
                        const RangeType& bval,
                        ResultType& resEn) {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      ResultType dummy(0.0);
      const Geometry& global = nit.intersectionGlobal();
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(quad.point(quadPoint)),
                    valuesEn_);
      return problem_.numericalFlux(nit, 0.0,
                                    global.global(quad.point(quadPoint)),
                                    valuesEn_.second, bval,
                                    resEn, dummy);
    }
    */

    template <class Entity, class DomainType, class ResultType>
    void source(Entity& en, const DomainType& x, ResultType& res) {
      evaluateLocal(en, x, valuesEn_);
      evaluateJacobianLocal(en, x);
      problem_.source(en, 0.0, x, valuesEn_.second, jacobians_, res);
      //CallerType::source(problem_, en, x, valuesEn_.second, jacobians_, res);
    }

    template <class Entity, class QuadratureType, class ResultType>
    void source(Entity& en, QuadratureType& quad, int quadPoint, 
                ResultType& res) 
    {
      evaluateQuad(en, quad, quadPoint, valuesEn_);
      evaluateJacobianQuad(en, quad, quadPoint);
      problem_.source(en, 0.0, quad.point(quadPoint), valuesEn_.second,
                      jacobians_, res);
      //CallerType::source(problem_, en, quad, quadPoint,
      //valuesEn_.second, jacobians_, res);
    }

  private:
    template <class Entity>
    void setter(Entity& en, LocalFunctionTupleType& tuple) {
      ForEachValuePair<DiscreteFunctionTupleType, 
        LocalFunctionTupleType> forEach(discreteFunctions_, tuple);
      LocalFunctionSetter<Entity> setter(en);
      forEach.apply(setter);
    }

    template <class Entity, class DomainType>
    void evaluateLocal(Entity& en, const DomainType& x, ValuePair& p) {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(p.first,
                                                        p.second);
      
      LocalFunctionEvaluateLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    template <class Entity, class QuadratureType>
    void evaluateQuad(Entity& en, QuadratureType& quad, int quadPoint, 
                      ValuePair& p) {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(p.first,
                                                        p.second);
      LocalFunctionEvaluateQuad<Entity, QuadratureType> eval(en, quad, quadPoint);
      forEach.apply(eval);
    }

    template <class Entity, class DomainType>
    void evaluateJacobianLocal(Entity& en, const DomainType& x)
    {
      ForEachValuePair<
        LocalFunctionTupleType,JacobianRangeTupleType> forEach(valuesEn_.first,
                                                               jacobians_);
      
      LocalFunctionEvaluateJacobianLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    template <class Entity, class Quadrature>
    void evaluateJacobianQuad(Entity& en, Quadrature& quad, int quadPoint) {
      ForEachValuePair<
        LocalFunctionTupleType,JacobianRangeTupleType> forEach(valuesEn_.first,
                                                               jacobians_);
      
      LocalFunctionEvaluateJacobianQuad<Entity, Quadrature> eval(en, quad, quadPoint);
      forEach.apply(eval);
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
    JacobianRangeTupleType jacobians_;
  };

}

#endif
