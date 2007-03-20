#ifndef DUNE_ELLIPTICDISCRETEMODELCALLER_HH
#define DUNE_ELLIPTICDISCRETEMODELCALLER_HH

#include <utility>
#include <memory>

#include <dune/common/fvector.hh> 

#include "callerutility.hh"
#include "ellipticdiscretemodel.hh"
#include "modelcaller.hh"

#include <dune/fem/misc/boundaryidentifier.hh>

namespace Dune {

  /**
   * @brief Wrapper class for all the template magic used to call the problem
   * methods.
   */
  template <class DiscreteModelImp, class ArgumentImp, class SelectorImp>
  class EllipticDiscreteModelCaller 
    : public DiscreteModelCaller<DiscreteModelImp,ArgumentImp,SelectorImp> 
  {
  public:
    typedef DiscreteModelCaller<DiscreteModelImp,ArgumentImp,SelectorImp> BaseType; 
    typedef DiscreteModelImp DiscreteModelType;
    typedef ArgumentImp TotalArgumentType;
    typedef typename SelectorImp::Base SelectorType;

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

    typedef BoundaryIdentifier BoundaryIdentifierType;
  public:
    EllipticDiscreteModelCaller(DiscreteModelType& problem) 
      : BaseType(problem) 
      , problem_(problem)
    {}

    // Ensure: entities set correctly before call
    template <class QuadratureType, class CoefficientType>
    void evaluateCoefficientFace(const IntersectionIterator& nit,
                         const QuadratureType& quadInner, 
                         const QuadratureType& quadOuter, 
                         const int quadPoint,
                         CoefficientType& coeffLeft, 
                         CoefficientType& coeffRight) 
    {
      // evaluate data functions 
      this->evaluateQuad( quadInner, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      this->evaluateQuad( quadOuter, quadPoint,
                         this->data_->localFunctionsNeigh(), this->valuesNeigh_);

      problem_.coefficientFace(nit, this->time_, 
                               quadInner.localPoint(quadPoint),
                               this->valuesEn_, 
                               this->valuesNeigh_,
                               coeffLeft,coeffRight);
    }

    // Ensure: entities set correctly before call
    template <class QuadratureType, class CoefficientType>
    void evaluateCoefficientBoundary(const IntersectionIterator& nit,
                         const QuadratureType& quadInner, 
                         const int quadPoint,
                         CoefficientType& coeff) 
    {
      this->evaluateQuad( quadInner, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      problem_.coefficient(this->data_->self(), 
                           this->time_, quadInner.point(quadPoint), 
                           this->valuesEn_, coeff);
    }
      
    template <class CoefficientType>
    void evaluateCoefficient(Entity& en, VolumeQuadratureType& quad, int quadPoint,
                        CoefficientType& coeff) 
    {
      evaluateQuad(quad, quadPoint, 
                   this->data_->localFunctionsSelf(), this->valuesEn_);
      problem_.coefficient(en, this->time_, quad.point(quadPoint), 
                           this->valuesEn_, coeff);
    }

    template <class IntersectionIterator>
    BoundaryIdentifierType
    boundaryValue(IntersectionIterator& nit,
                         const FaceQuadratureType& quad,
                         const int quadPoint,
                         RangeType& boundaryValue) 
    {
      typedef typename IntersectionIterator::LocalGeometry Geometry;

      this->evaluateQuad( quad, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      return problem_.boundaryValue(nit, this->time_, quad.localPoint(quadPoint),
                                    this->valuesEn_, boundaryValue);
    }

    void rightHandSide(Entity& en, VolumeQuadratureType& quad, int quadPoint, 
                       RangeType& res) 
    {
      this->evaluateQuad( quad, quadPoint,
                         this->data_->localFunctionsSelf(), this->valuesEn_);
      this->evaluateJacobianQuad(en, quad, quadPoint);
      problem_.rightHandSide(en, this->time_, quad.point(quadPoint), 
                             this->valuesEn_,
                             this->jacobians_, res);
    }

  private:
    // our problem 
    DiscreteModelType& problem_;
  }; // end EllipticDiscreteModelCaller 

} // end namespace Dune 
#endif
