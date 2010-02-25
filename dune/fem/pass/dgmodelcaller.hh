#ifndef DUNE_DGDISCRETEMODELCALLER_HH
#define DUNE_DGDISCRETEMODELCALLER_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include "callerutility.hh"
#include "dgdiscretemodel.hh"
#include "modelcallerdefault.hh"

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class DGDiscreteModelCaller 
    : public DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp >
  {
  public:
    typedef DiscreteModelCallerDefault< DiscreteModelImp, ArgumentImp, SelectorImp > BaseType;
    typedef typename BaseType::DataStorage DataStorage;
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
    typedef typename IntersectionIterator::Intersection IntersectionType;
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
    DGDiscreteModelCaller( DiscreteModelType& problem ) :
      BaseType( ), 
      problem_( problem )
    {}

    void setEntity ( const Entity &entity ) 
    {
      BaseType::setEntity( entity );
      problem_.setEntity( entity );
    }

    void setNeighbor( const Entity &neighbor )
    {
      BaseType::setNeighbor( neighbor );
      problem_.setNeighbor( neighbor );
    }

    // Here, the interface of the problem is replicated and the Caller
    // is used to do the actual work
    void analyticalFlux( const Entity& en, 
                         const DomainType& x, 
                         JacobianRangeType& res) 
    {
      evaluateLocal(x, data_->localFunctionsSelf(), valuesEn_);
      problem_.analyticalFlux(en, time_, x, valuesEn_, res);
    }
    
    void analyticalFlux ( const Entity &entity,
                          const VolumeQuadratureType &quad,  const int quadPoint,
                          JacobianRangeType &res )
    {
      evaluateQuad( quad, quadPoint, data_->localFunctionsSelf(), valuesEn_ );
      problem_.analyticalFlux( entity, time_, quad.point( quadPoint ), valuesEn_, res );
    }

    void analyticalFluxAndSource( const Entity &entity,
                                  const VolumeQuadratureType& quad, 
                                  const int quadPoint,
                                  JacobianRangeType& fluxRes, 
                                  RangeType& sourceRes ) 
    {
      // evaluate local functions 
      evaluateQuad(quad, quadPoint, data_->localFunctionsSelf(),valuesEn_);
      // evaluate local gradients 
      evaluateJacobianQuad( quad, quadPoint, data_->localFunctionsSelf(), jacobians_ );
      
      problem_.analyticalFlux( entity, time_, quad.point( quadPoint ), valuesEn_, fluxRes );
      problem_.source( entity, time_, quad.point( quadPoint ), valuesEn_, jacobians_, sourceRes );
    }

    template< class IntersectionIterator, class FaceDomainType >
    double numericalFlux( const IntersectionIterator& nit,
                          const FaceDomainType& x,
                          RangeType& resEn,
                          RangeType& resNeigh ) 
    {
      typedef typename IntersectionIterator::Intersection::LocalGeometry Geometry;

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
    template< class IntersectionIterator, class QuadratureType >
    double numericalFlux( const IntersectionIterator& nit,
                          const QuadratureType& quadInner, 
                          const QuadratureType& quadOuter, 
                          const int quadPoint,
                          RangeType& resEn,
                          RangeType& resNeigh )
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
    template< class IntersectionIterator, class QuadratureType >
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
      typedef typename IntersectionIterator::Intersection::LocalGeometry Geometry;
      const Geometry& selfLocal = nit.intersectionSelfLocal();
      evaluateLocal(selfLocal.global(x),
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, x,
                                   valuesEn_,
                                   boundaryFlux);
    }

    //template <class IntersectionIterator>
    double boundaryFlux(const IntersectionType& nit,
                        const FaceQuadratureType& quad,
                        const int quadPoint,
                        RangeType& boundaryFlux) 
    {
      typedef typename IntersectionIterator::Intersection::LocalGeometry Geometry;

      evaluateQuad( quad, quadPoint,
                    data_->localFunctionsSelf(), valuesEn_);
      return problem_.boundaryFlux(nit, time_, quad.localPoint(quadPoint),
                                   valuesEn_, boundaryFlux);
    }

    void source( const Entity& en, 
                 const DomainType& x, 
                 RangeType& res ) 
    {
      evaluateLocal( x, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianLocal(en, x);
      problem_.source(en, time_, x, valuesEn_, jacobians_, res);
    }

    void source( const Entity& en, 
                 const VolumeQuadratureType& quad, 
                 const int quadPoint, 
                 RangeType& res ) 
    {
      evaluateQuad( quad, quadPoint, data_->localFunctionsSelf(), valuesEn_);
      evaluateJacobianQuad( quad, quadPoint, data_->localFunctionsSelf(), jacobians_);
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
    
  private:
    DGDiscreteModelCaller(const DGDiscreteModelCaller&);
    DGDiscreteModelCaller& operator=(const DGDiscreteModelCaller&);

  protected:
    DiscreteModelType& problem_;

  protected:  
    using BaseType::data_;
    using BaseType::valuesEn_;
    using BaseType::valuesNeigh_;
    using BaseType::jacobians_;
    using BaseType::time_;
  };

}

#endif
