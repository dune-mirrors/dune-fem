#ifndef DUNE_DISCRETEMODELCALLERDEFAULT_HH
#define DUNE_DISCRETEMODELCALLERDEFAULT_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include "callerutility.hh"

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class DiscreteModelCallerDefault
  {
  public:
    typedef DiscreteModelImp DiscreteModelType;
    typedef ArgumentImp TotalArgumentType;
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
    DiscreteModelCallerDefault(  ) :
      data_( 0 ),
      valuesEn_( RangeCreator::apply() ),
      valuesNeigh_( RangeCreator::apply() ),
      jacobians_( JacobianCreator::apply() ),
      time_( 0.0 )
    {}

    void setArgument(TotalArgumentType& arg) 
    {
      data_.reset(new DataStorage(arg));
    }

    void setEntity(Entity& en) 
    {
      data_->setSelf(en);
    }

    void setNeighbor(Entity& nb) 
    {
      data_->setNeighbor(nb);
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

  protected:
    void setter ( Entity &entity, LocalFunctionTupleType &localFunctionTuple )
    {
      ForEachValue< LocalFunctionTupleType > forEach( localFunctionTuple );
      LocalFunctionSetter< Entity > setter( entity );
      forEach.apply( setter );
    }

#if 0
    void evaluateLocal(const DomainType& x, 
                       LocalFunctionTupleType& lfs, RangeTupleType& ranges) 
    {
      ForEachValuePair<
        LocalFunctionTupleType, RangeTupleType> forEach(lfs,
                                                        ranges);
      
      LocalFunctionEvaluateLocal<DomainType> eval(x);
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

#endif

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
    DiscreteModelCallerDefault(const DiscreteModelCallerDefault&);
    DiscreteModelCallerDefault& operator=(const DiscreteModelCallerDefault&);

  protected:
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

  protected:  
    std::auto_ptr<DataStorage> data_;

    RangeTupleType valuesEn_;
    RangeTupleType valuesNeigh_;
    JacobianRangeTupleType jacobians_;

    double time_;
  };

}

#endif
