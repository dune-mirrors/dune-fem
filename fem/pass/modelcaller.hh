#ifndef DUNE_DISCRETEMODELCALLER_HH
#define DUNE_DISCRETEMODELCALLER_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include "callerutility.hh"
#include "discretemodel.hh"

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class DiscreteModelCaller {
  public:
    typedef DiscreteModelImp DiscreteModelType;
    typedef ArgumentImp TotalArgumentType;
    typedef typename SelectorImp::Base SelectorType;

    typedef typename DiscreteModelType::Traits Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridType GridType;
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity Entity;

    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;

    typedef Filter<TotalArgumentType, SelectorType> FilterType;
    typedef typename FilterType::ResultType DiscreteFunctionTupleType;
    typedef LocalFunctionCreator<DiscreteFunctionTupleType> LFCreator;
    typedef typename LFCreator::ResultType LocalFunctionTupleType;
    typedef Creator<
      RangeTypeEvaluator, LocalFunctionTupleType> RangeCreator;
    typedef typename RangeCreator::ResultType RangeTupleType;
    typedef Creator<
      JacobianRangeTypeEvaluator, LocalFunctionTupleType> JacobianCreator;
    typedef typename JacobianCreator::ResultType JacobianRangeTupleType;
    //typedef typename Caller<RangeTupleType> CallerType;

  public:
    DiscreteModelCaller(DiscreteModelType& problem) :
      problem_(problem),
      data_(0),
      valuesEn_(RangeCreator::apply()),
      valuesNeigh_(RangeCreator::apply()),
      jacobians_(JacobianCreator::apply()),
      time_(0.0)
    {}

    void setArgument(TotalArgumentType& arg) 
    {
      data_.reset(new DataStorage(arg));
    }

    void setEntity(Entity& en) 
    {
      data_->setSelf(en);
    }

    void setNeighbor(Entity& en) 
    {
      data_->setNeighbor(en);
    }

    void setTime(double time) {
      time_ = time;
    }

    void finalize() {
      data_.reset(0);
    }

    template <class QuadratureType>
    void setQuad(Entity& en,QuadratureType& quad) 
    {
      ForEachValue<DiscreteFunctionTupleType> forEach(data_->discreteFunctions());
      DiscreteFunctionSetQuad<Entity,QuadratureType> eval(en,quad);
      forEach.apply(eval);
    }

    // Here, the interface of the problem is replicated and the Caller
    // is used to do the actual work
    void analyticalFlux(Entity& en, const DomainType& x, 
                        JacobianRangeType& res) 
    {
      evaluateLocal(en, x, data_->localFunctionsSelf(), valuesEn_);
      problem_.analyticalFlux(en, time_, x, valuesEn_, res);
      //CallerType::analyticalFlux(problem_, en, x, valuesEn_.second, res);
    }
    
    void analyticalFlux(Entity& en, VolumeQuadratureType& quad, int quadPoint,
                        JacobianRangeType& res) 
    {
      evaluateQuad(en, quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      problem_.analyticalFlux(en, time_, quad.point(quadPoint), 
                              valuesEn_, res);
      //CallerType::analyticalFlux(problem_, en, quad, quadPoint, valuesEn_.second, res);
    }

    template <class FaceDomainType>
    double numericalFlux(IntersectionIterator& nit, const FaceDomainType& x,
                         RangeType& resEn, RangeType& resNeigh) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;

      const Geometry& selfLocal = nit.intersectionSelfLocal();
      const Geometry& neighLocal = nit.intersectionNeighborLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(x),
                    data_->localFunctionsSelf(), valuesEn_);
      evaluateLocal(*nit.outside(), neighLocal.global(x),
                    data_->localFunctionsNeigh(), valuesNeigh_);
      return problem_.numericalFlux(nit, time_, x,
                                    valuesEn_, valuesNeigh_,
                                    resEn, resNeigh);
      //CallerType::numericalFlux(problem_, nit, x, 
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }

    // Ensure: entities set correctly before call
    double numericalFlux(IntersectionIterator& nit,
                         FaceQuadratureType& quadInner, 
                         FaceQuadratureType& quadOuter, 
                         int quadPoint,
                         RangeType& resEn, RangeType& resNeigh)
    {
      evaluateQuad(*nit.inside(), quadInner, quadPoint,
                   data_->localFunctionsSelf(), valuesEn_);
      evaluateQuad(*nit.outside(), quadOuter, quadPoint,
                   data_->localFunctionsNeigh(), valuesNeigh_);
      return problem_.numericalFlux(nit, time_, 
                                    quadInner.localPoint(quadPoint),
                                    valuesEn_, 
                                    valuesNeigh_,
                                    resEn, resNeigh);
      //CallerType::numericalFlux(problem_, nit, quad, quadPoint,
      //                      valuesEn_.second, valuesNeigh_.second,
      //                      res_en, res_neigh);
    }


    template <class IntersectionIterator, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& nit,
                        const FaceDomainType& x,
                        RangeType& boundaryFlux) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(x),
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, x,
                                   valuesEn_,
                                   boundaryFlux);
    }

    template <class IntersectionIterator>
    double boundaryFlux(IntersectionIterator& nit,
                        FaceQuadratureType& quad,
                        int quadPoint,
                        RangeType& boundaryFlux) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;

      evaluateQuad(*nit.inside(), quad, quadPoint,
                   data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, quad.localPoint(quadPoint),
                                   valuesEn_, boundaryFlux);
    }

    void source(Entity& en, const DomainType& x, RangeType& res) 
    {
      evaluateLocal(en, x, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianLocal(en, x);
      problem_.source(en, time_, x, valuesEn_, jacobians_, res);
      //CallerType::source(problem_, en, x, valuesEn_.second, jacobians_, res);
    }

    void source(Entity& en, VolumeQuadratureType& quad, int quadPoint, 
                RangeType& res) 
    {
      // evaluateQuad(en, quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      evaluateJacobianQuad(en, quad, quadPoint);
      problem_.source(en, time_, quad.point(quadPoint), valuesEn_,
                      jacobians_, res);
      //CallerType::source(problem_, en, quad, quadPoint,
      //valuesEn_.second, jacobians_, res);
    }

  private:
    void setter(Entity& en, LocalFunctionTupleType& tuple) 
    {
      ForEachValue<LocalFunctionTupleType> forEach(tuple);
      LocalFunctionSetter<Entity> setter(en);
      forEach.apply(setter);
    }

    void evaluateLocal(Entity& en, const DomainType& x, 
                       LocalFunctionTupleType& lfs, RangeTupleType& ranges) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(lfs,
                                                        ranges);
      
      LocalFunctionEvaluateLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    template <class QuadratureType>
    inline
    void evaluateQuad(Entity& en, QuadratureType& quad, int quadPoint, 
                      LocalFunctionTupleType& lfs, RangeTupleType& ranges) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(lfs,
                                                        ranges);
      LocalFunctionEvaluateQuad<
        Entity, QuadratureType> eval(en, quad, quadPoint);
      
      forEach.apply(eval);
    }

    void evaluateJacobianLocal(Entity& en, const DomainType& x)
    {
      ForEachValuePair<
       LocalFunctionTupleType,
        JacobianRangeTupleType> forEach(data_->localFunctionsSelf(),
                                        jacobians_);
      
      LocalFunctionEvaluateJacobianLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    template <class QuadratureType>
    void evaluateJacobianQuad(Entity& en, QuadratureType& quad, 
                              int quadPoint) 
    {
      ForEachValuePair<
       LocalFunctionTupleType,
        JacobianRangeTupleType> forEach(data_->localFunctionsSelf(),
                                        jacobians_);
      LocalFunctionEvaluateJacobianQuad<
        Entity, QuadratureType> eval(en, quad, quadPoint);
      
      forEach.apply(eval);
    }

  private:
    DiscreteModelCaller(const DiscreteModelCaller&);
    DiscreteModelCaller& operator=(const DiscreteModelCaller&);

  private:
    class DataStorage 
    {
    public:
      DataStorage(TotalArgumentType& arg) :
        discreteFunctions_(FilterType::apply(arg)),
        localFunctionsSelf_(LFCreator::apply(discreteFunctions_)),
        localFunctionsNeigh_(LFCreator::apply(discreteFunctions_))
      {}

      LocalFunctionTupleType& localFunctionsSelf() {
        return localFunctionsSelf_;
      }

      LocalFunctionTupleType& localFunctionsNeigh() {
        return localFunctionsNeigh_;
      }

      DiscreteFunctionTupleType& discreteFunctions() {
	return discreteFunctions_;
      }

      void setSelf(Entity& en) {
        setter(en, localFunctionsSelf_);
      }

      void setNeighbor(Entity& en) {
        setter(en, localFunctionsNeigh_);
      }

    private:
      void setter(Entity& en, LocalFunctionTupleType& tuple) {
        ForEachValue<LocalFunctionTupleType> forEach(tuple);
        LocalFunctionSetter<Entity> setter(en);
        forEach.apply(setter);
      }

    private:
      DiscreteFunctionTupleType discreteFunctions_;
      LocalFunctionTupleType localFunctionsSelf_;
      LocalFunctionTupleType localFunctionsNeigh_;
    };

  private:
    DiscreteModelType& problem_;
    
    std::auto_ptr<DataStorage> data_;

    RangeTupleType valuesEn_;
    RangeTupleType valuesNeigh_;
    JacobianRangeTupleType jacobians_;

    double time_;
  };

}

#endif
