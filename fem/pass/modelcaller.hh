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
    //typedef typename SelectorImp::Base SelectorType;
    typedef SelectorImp SelectorType;

    typedef typename DiscreteModelType::Traits Traits;
    typedef typename Traits::DomainType DomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename GridPartType :: GridType GridType;
    typedef typename GridPartType::IntersectionIteratorType IntersectionIterator;
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

    //! type of mass matrix factor (see discretemodel.hh)
    typedef typename DiscreteModelType :: MassFactorType MassFactorType;
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
      problem_.setEntity(en);
    }

    void setNeighbor(Entity& nb) 
    {
      data_->setNeighbor(nb);
      problem_.setNeighbor(nb);
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

    template <class QuadratureType>
    void setQuadSelf(QuadratureType& quad) 
    {
      ForEachValue<LocalFunctionTupleType> forEach(data_->localFunctionsSelf());
      LocalDiscreteFunctionSetQuad<QuadratureType> eval(quad);
      forEach.apply(eval);
    }

    template <class QuadratureType>
    void setQuadNeigh(QuadratureType& quad) 
    {
      ForEachValue<LocalFunctionTupleType> forEach(data_->localFunctionsNeigh());
      LocalDiscreteFunctionSetQuad<QuadratureType> eval(quad);
      forEach.apply(eval);
    }

    // Here, the interface of the problem is replicated and the Caller
    // is used to do the actual work
    void analyticalFlux(Entity& en, const DomainType& x, 
                        JacobianRangeType& res) 
    {
      evaluateLocal(x, data_->localFunctionsSelf(), valuesEn_);
      problem_.analyticalFlux(en, time_, x, valuesEn_, res);
    }
    
    void analyticalFlux(Entity& en, VolumeQuadratureType& quad, int quadPoint,
                        JacobianRangeType& res) 
    {
      evaluateQuad(quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      problem_.analyticalFlux(en, time_, quad.point(quadPoint), 
                              valuesEn_, res);
    }

    void analyticalFluxAndSource(Entity& en, VolumeQuadratureType& quad, int quadPoint,
                                 JacobianRangeType& fluxRes, RangeType& sourceRes) 
    {
      // evaluate local functions 
      evaluateQuad(quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      // evaluate local gradients 
      evaluateJacobianQuad(en, quad, quadPoint);
      
      problem_.analyticalFlux(en, time_, quad.point(quadPoint), 
                              valuesEn_, fluxRes);

      problem_.source(en, time_, quad.point(quadPoint), valuesEn_,
                      jacobians_, sourceRes);
    }

    template <class FaceDomainType>
    double numericalFlux(IntersectionIterator& nit, const FaceDomainType& x,
                         RangeType& resEn, RangeType& resNeigh) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;

      const Geometry& selfLocal = nit.intersectionSelfLocal();
      const Geometry& neighLocal = nit.intersectionNeighborLocal();
      evaluateLocal( selfLocal.global(x),
                    data_->localFunctionsSelf(), valuesEn_);
      evaluateLocal( neighLocal.global(x),
                    data_->localFunctionsNeigh(), valuesNeigh_);
      return problem_.numericalFlux(nit, time_, x,
                                    valuesEn_, valuesNeigh_,
                                    resEn, resNeigh);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType>
    double numericalFlux(IntersectionIterator& nit,
                         QuadratureType& quadInner, 
                         QuadratureType& quadOuter, 
                         int quadPoint,
                         RangeType& resEn, RangeType& resNeigh)
    {
      evaluateQuad(quadInner, quadPoint,
                   data_->localFunctionsSelf(), valuesEn_);
      evaluateQuad(quadOuter, quadPoint,
                   data_->localFunctionsNeigh(), valuesNeigh_);
      return problem_.numericalFlux(nit, time_, 
                                    quadInner.localPoint(quadPoint),
                                    valuesEn_, 
                                    valuesNeigh_,
                                    resEn, resNeigh);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType>
    double numericalFlux(const IntersectionIterator& nit,
                         const QuadratureType& quadInner, 
                         const QuadratureType& quadOuter, 
                         const int quadPoint,
                         RangeType & sigmaPartLeft, 
                         RangeType & sigmaPartRight,
                         RangeType& uPartLeft, 
                         RangeType& uPartRight)
    {
      // evaluate data functions 
      evaluateQuad( quadInner, quadPoint,
                    data_->localFunctionsSelf(), valuesEn_);
      evaluateQuad( quadOuter, quadPoint,
                    data_->localFunctionsNeigh(), valuesNeigh_);

      return problem_.numericalFlux(nit, time_, 
                                    quadInner.localPoint(quadPoint),
                                    valuesEn_, 
                                    valuesNeigh_,
                                    sigmaPartLeft,sigmaPartRight,
                                    uPartLeft, uPartRight);
    }

    template <class IntersectionIterator, class FaceDomainType>
    double boundaryFlux(IntersectionIterator& nit,
                        const FaceDomainType& x,
                        RangeType& boundaryFlux) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(selfLocal.global(x),
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, x,
                                   valuesEn_,
                                   boundaryFlux);
    }

    template <class IntersectionIterator>
    double boundaryFlux(IntersectionIterator& nit,
                        const FaceQuadratureType& quad,
                        const int quadPoint,
                        RangeType& boundaryFlux) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;

      evaluateQuad( quad, quadPoint,
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, quad.localPoint(quadPoint),
                                   valuesEn_, boundaryFlux);
    }

    void source(Entity& en, const DomainType& x, RangeType& res) 
    {
      evaluateLocal( x, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianLocal(en, x);
      problem_.source(en, time_, x, valuesEn_, jacobians_, res);
    }

    void source(Entity& en, VolumeQuadratureType& quad, int quadPoint, 
                RangeType& res) 
    {
      evaluateQuad( quad, quadPoint, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianQuad(en, quad, quadPoint);
      problem_.source(en, time_, quad.point(quadPoint), valuesEn_,
                      jacobians_, res);
    }

    //! return true when a mass matrix has to be build  
    bool hasMass () const  { return problem_.hasMass(); }

    //! evaluate mass matrix factor 
    void mass(const Entity& en,
              const VolumeQuadratureType& quad,
              const int quadPoint,
              MassFactorType& m)
    {
      // evaluate local functions 
      evaluateQuad( quad, quadPoint,
                    data_->localFunctionsSelf(), valuesEn_);
      // call problem implementation 
      problem_.mass(en, time_, quad.point(quadPoint),
                    valuesEn_,
                    m);
    }
    
  protected:
    void setter ( Entity &entity, LocalFunctionTupleType &localFunctionTuple )
    {
      ForEachValue< LocalFunctionTupleType > forEach( localFunctionTuple );
      LocalFunctionSetter< Entity > setter( entity );
      forEach.apply( setter );
    }

    void evaluateLocal(const DomainType& x, 
                       LocalFunctionTupleType& lfs, RangeTupleType& ranges) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(lfs,
                                                        ranges);
      
      LocalFunctionEvaluateLocal<DomainType> eval(x);
      forEach.apply(eval);
    }

    template <class QuadratureType>
    inline
    void evaluateQuad(const QuadratureType& quad, 
                      const int quadPoint, 
                      LocalFunctionTupleType& lfs, 
                      RangeTupleType& ranges) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(lfs,
                                                        ranges);
      LocalFunctionEvaluateQuad<QuadratureType> eval(quad, quadPoint);
      
      forEach.apply(eval);
    }

    void evaluateJacobianLocal(Entity& en, const DomainType& x)
    {
      ForEachValuePair<
       LocalFunctionTupleType,
        JacobianRangeTupleType> forEach(data_->localFunctionsSelf(),
                                        jacobians_);
      
      LocalFunctionEvaluateJacobianLocal<DomainType> eval(x);
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
        QuadratureType> eval(quad, quadPoint);
      
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
        localFunctionsNeigh_(LFCreator::apply(discreteFunctions_)),
      	self_(0),
      	neighbor_(0)
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

      Entity& self() {
        assert( self_ );
        return *self_;
      }

      Entity& neighbor() {
        assert( neighbor_ );
        return *neighbor_;
      }

      void setSelf(Entity& en) {
      	self_ = &en;
        setter(en, localFunctionsSelf_);
      }

      void setNeighbor(Entity& en) {
      	neighbor_=&en;
        setter(en, localFunctionsNeigh_);
      }

    private:
      void setter( Entity &entity, LocalFunctionTupleType &localFunctionTuple )
      {
        ForEachValue< LocalFunctionTupleType > forEach( localFunctionTuple );
        LocalFunctionSetter< Entity > setter( entity );
        forEach.apply( setter );
      }

    private:
      DiscreteFunctionTupleType discreteFunctions_;
      LocalFunctionTupleType localFunctionsSelf_;
      LocalFunctionTupleType localFunctionsNeigh_;

      Entity* self_;
      Entity* neighbor_;
    };

  private:
    DiscreteModelType& problem_;

  protected:  
    std::auto_ptr<DataStorage> data_;

    RangeTupleType valuesEn_;
    RangeTupleType valuesNeigh_;
    JacobianRangeTupleType jacobians_;

    // Entity* self_;
    // Entity* neighbor_;
    
    double time_;
  };

}

#endif
