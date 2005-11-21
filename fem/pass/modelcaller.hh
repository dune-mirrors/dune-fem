#ifndef DUNE_DISCRETEMODELCALLER_HH
#define DUNE_DISCRETEMODELCALLER_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include "callerutility.hh"
#include "discretemodel.hh"

namespace Dune {

  //! \todo Provide shortcuts for evaluation on quadrature points once the caching is in place

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class DiscreteModelCaller {
  public:
    typedef DiscreteModelImp DiscreteModelType;
    typedef ArgumentImp TotalArgumentType;
    //typedef SelectorImp SelectorType;
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
      assert(data_.get());
      data_->setSelf(en);
    }

    void setNeighbor(Entity& en) 
    {
      assert(data_.get());
      data_->setNeighbor(en);
    }

    void setTime(double time) {
      time_ = time;
    }

    void finalize() {
      data_.reset(0);
    }

    // Here, the interface of the problem is replicated and the Caller
    // is used to do the actual work
    void analyticalFlux(Entity& en, const DomainType& x, 
                        JacobianRangeType& res) 
    {
      assert(data_.get());

      evaluateLocal(en, x, data_->localFunctionsSelf(), valuesEn_);
      problem_.analyticalFlux(en, time_, x, valuesEn_, res);
      //CallerType::analyticalFlux(problem_, en, x, valuesEn_.second, res);
    }
    
    void analyticalFlux(Entity& en, VolumeQuadratureType& quad, int quadPoint,
                        JacobianRangeType& res) 
    {
      // * temporary
      /*
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(*valuesEn_.first,
                                                        valuesEn_.second);
      LocalFunctionEvaluateQuad<
        Entity, VolumeQuadratureType> eval(en, 
                                           quad, 
                                           quadPoint);
      forEach.apply(eval);
      */
      assert(data_.get());

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

      assert(data_.get());

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
      //                      valuesEn_, valuesNeigh_,
      //                      res_en, res_neigh);
    }

    // Ensure: entities set correctly before call
    double numericalFlux(IntersectionIterator& nit,
                         FaceQuadratureType& quad, int quadPoint,
                         RangeType& resEn, RangeType& resNeigh)
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
    
      assert(data_.get());

      const Geometry& selfLocal = nit.intersectionSelfLocal();
      const Geometry& neighLocal = nit.intersectionNeighborLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(quad.point(quadPoint)),
                    data_->localFunctionsSelf(), valuesEn_);
      evaluateLocal(*nit.outside(), neighLocal.global(quad.point(quadPoint)),
                    data_->localFunctionsNeigh(), valuesNeigh_);
      return problem_.numericalFlux(nit, time_, 
                                    quad.point(quadPoint),
                                    valuesEn_, valuesNeigh_,
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

      assert(data_.get());

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

      assert(data_.get());

      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(*nit.inside(), selfLocal.global(quad.point(quadPoint)),
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, quad.point(quadPoint),
                                   valuesEn_,
                                   boundaryFlux);
    }

    void source(Entity& en, const DomainType& x, RangeType& res) 
    {
      assert(data_.get());

      evaluateLocal(en, x, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianLocal(en, x);
      problem_.source(en, time_, x, valuesEn_, jacobians_, res);
      //CallerType::source(problem_, en, x, valuesEn_, jacobians_, res);
    }

    void source(Entity& en, VolumeQuadratureType& quad, int quadPoint, 
                RangeType& res) 
    {
      assert(data_.get());

      evaluateQuad(en, quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      evaluateJacobianQuad(en, quad, quadPoint);
      problem_.source(en, time_, quad.point(quadPoint), valuesEn_,
                      jacobians_, res);
      //CallerType::source(problem_, en, quad, quadPoint,
      //valuesEn_.second, jacobians_, res);
    }

  private:
    void evaluateLocal(Entity& en, const DomainType& x, 
                       LocalFunctionTupleType& l, RangeTupleType& r) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(l, r);
      
      LocalFunctionEvaluateLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    void evaluateQuad(Entity& en, VolumeQuadratureType& quad, int quadPoint, 
                      LocalFunctionTupleType& l, RangeTupleType& r) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(l, r);

      LocalFunctionEvaluateQuad<
        Entity, VolumeQuadratureType> eval(en, quad, quadPoint);
      
      forEach.apply(eval);
    }

    void evaluateJacobianLocal(Entity& en, const DomainType& x)
    {
      assert(data_.get());

      ForEachValuePair<LocalFunctionTupleType, JacobianRangeTupleType> 
        forEach(data_->localFunctions(), jacobians_);
      LocalFunctionEvaluateJacobianLocal<Entity, DomainType> eval(en, x);
      forEach.apply(eval);
    }

    void evaluateJacobianQuad(Entity& en, VolumeQuadratureType& quad, 
                              int quadPoint) 
    {
      assert(data_.get());

      ForEachValuePair<LocalFunctionTupleType, JacobianRangeTupleType> 
        forEach(data_->localFunctionsSelf(), jacobians_);
      LocalFunctionEvaluateJacobianQuad<
        Entity, VolumeQuadratureType> eval(en, quad, quadPoint);
      
      forEach.apply(eval);
    }

  private:
    DiscreteModelCaller(const DiscreteModelCaller&);
    DiscreteModelCaller& operator=(const DiscreteModelCaller&);

  private:
    /*
    DiscreteModelType& problem_;
    TotalArgumentType* arg_;

    DiscreteFunctionTupleType discreteFunctions_;
    ValuePair valuesEn_;
    ValuePair valuesNeigh_;
    JacobianRangeTupleType jacobians_;

    double time_;
    */

    class DataStorage {
    public:
      DataStorage(TotalArgumentType& arg) :
        discreteFunctions_(FilterType::apply(arg)),
        localFunctionsEn_(LFCreator::apply(discreteFunctions_)),
        localFunctionsNeigh_(LFCreator::apply(discreteFunctions_))
      {}

      DiscreteFunctionTupleType& discreteFunctions() 
      {
        return discreteFunctions_;
      }

      LocalFunctionTupleType& localFunctionsSelf()
      {
        return localFunctionsEn_;
      }

      LocalFunctionTupleType& localFunctionsNeigh()
      {
        return localFunctionsEn_;
      }

      void setSelf(Entity& en)
      {
        setter(en, localFunctionsEn_);
      }

      void setNeighbor(Entity& en)
      {
        setter(en, localFunctionsNeigh_);
      }

    private:
      void setter(Entity& en, LocalFunctionTupleType& tuple) 
      {
        ForEachValuePair<DiscreteFunctionTupleType, 
          LocalFunctionTupleType> forEach(discreteFunctions_, tuple);
        LocalFunctionSetter<Entity> setter(en);
        forEach.apply(setter);
      }

    private:
      DiscreteFunctionTupleType discreteFunctions_;
      LocalFunctionTupleType localFunctionsEn_;
      LocalFunctionTupleType localFunctionsNeigh_;
    };

  private:
    //- Discrete model
    DiscreteModelType& problem_;
    
    //- DataStorage
    std::auto_ptr<DataStorage> data_;

    //- Data
    RangeTupleType valuesEn_;
    RangeTupleType valuesNeigh_;
    JacobianRangeTupleType jacobians_;

    //- Time
    double time_;
  };

}

#endif
