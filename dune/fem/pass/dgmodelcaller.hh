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
#ifndef NDEBUG
      , enQuadId_( (size_t) -1 ) 
      , faceQuadId_( (size_t) -1 ) 
#endif
    {
      // we need at least size 1 
      valuesEnVec_.push_back( RangeTupleType( RangeCreator::apply() ) );
      valuesNbVec_.push_back( RangeTupleType( RangeCreator::apply() ) );

      valuesVec_.push_back( RangeTupleType( RangeCreator::apply() ) );
      jacobianVec_.push_back( JacobianRangeTupleType( JacobianCreator::apply() )  );
    }

    void setEntity ( const Entity &entity ) 
    {
      BaseType::setEntity( entity );
      problem_.setEntity( entity );
    }

    template <class QuadratureType>
    void setEntity ( const Entity &entity, 
                     const QuadratureType& quad)
    {
      setEntity( entity );

      // evaluate all local functions for whole quadrature 
      evaluateQuadrature( quad, data_->localFunctionsSelf(), valuesVec_ );

      // only when we have a source term we need the jacobians 
      if( problem_.hasSource () ) 
      {
        // evaluate all local functions for whole quadrature 
        evaluateQuadrature( quad, 
            data_->localFunctionsSelf(), jacobianVec_ );
      }

#ifndef NDEBUG 
      enQuadId_ = quad.id();
#endif
    }

    void setNeighbor( const Entity &neighbor )
    {
      BaseType::setNeighbor( neighbor );
      problem_.setNeighbor( neighbor );
    }

    template< class QuadratureType >
    void setNeighbor( const Entity& neighbor,
                      const QuadratureType& quadInner,
                      const QuadratureType& quadOuter)
    {
      setNeighbor( neighbor );

      // evaluate all local functions for whole quadrature 
      evaluateQuadrature( quadInner, data_->localFunctionsSelf(), valuesEnVec_ );

      // evaluate all local functions for whole quadrature 
      evaluateQuadrature( quadOuter, data_->localFunctionsNeigh(), valuesNbVec_ );

#ifndef NDEBUG 
      faceQuadId_ = quadInner.id();
#endif
    }
    
    template< class QuadratureType >
    void setBoundary( const Entity& entity,
                      const QuadratureType& quad)
    {
      // evaluate all local functions for whole quadrature 
      evaluateQuadrature( quad, data_->localFunctionsSelf(), valuesEnVec_ );

#ifndef NDEBUG 
      faceQuadId_ = quad.id();
#endif
    }

    void analyticalFlux ( const Entity &entity,
                          const VolumeQuadratureType &quad,  
                          const int quadPoint,
                          JacobianRangeType &res )
    {
      // make sure we git the right quadrature 
      assert( enQuadId_ == quad.id() );
      assert( (int) valuesVec_.size() > quadPoint );

      problem_.analyticalFlux( entity, time_, quad.point( quadPoint ), valuesVec_[ quadPoint ] , res );
    }

    void analyticalFluxAndSource( const Entity &entity,
                                  const VolumeQuadratureType& quad, 
                                  const int quadPoint,
                                  JacobianRangeType& fluxRes, 
                                  RangeType& sourceRes ) 
    {
      assert( enQuadId_ == quad.id() );
      assert( (int) valuesVec_.size() > quadPoint );
      assert( (int) jacobianVec_.size() > quadPoint );

      problem_.analyticalFlux( entity, time_, quad.point( quadPoint ), 
          valuesVec_[ quadPoint ], fluxRes );
      problem_.source( entity, time_, quad.point( quadPoint ), 
          valuesVec_[ quadPoint ], jacobianVec_[ quadPoint ], sourceRes );
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
      assert( (int) valuesEnVec_.size() > quadPoint );
      // make sure we git the right quadrature 
      assert( faceQuadId_ == quadInner.id() );

      return problem_.numericalFlux(nit, time_, 
                                    quadInner.localPoint(quadPoint),
                                    valuesEnVec_[ quadPoint ], 
                                    valuesNbVec_[ quadPoint ],
                                    resEn, resNeigh);
    }


#if 0
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
      abort();
      // evaluate data functions 
      evaluateQuad( quadInner, quadPoint,
                    data_->localFunctionsSelf(), valuesEn_);
      evaluateQuad( quadOuter, quadPoint,
                    data_->localFunctionsNeigh(), valuesNb_);

      return problem_.numericalFlux(nit, time_, 
                                    quadInner.localPoint(quadPoint),
                                    valuesEn_, 
                                    valuesNb_,
                                    sigmaPartLeft,sigmaPartRight,
                                    uPartLeft, uPartRight);
    }
#endif

    double boundaryFlux(const IntersectionType& nit,
                        const FaceQuadratureType& quad,
                        const int quadPoint,
                        RangeType& boundaryFlux) 
    {
      assert( faceQuadId_ == quad.id() );
      assert( (int) valuesEnVec_.size() > quadPoint );

      return problem_.boundaryFlux(nit, time_, quad.localPoint(quadPoint),
                                   valuesEnVec_[ quadPoint ], boundaryFlux);
    }

    //! return true when a mass matrix has to be build  
    bool hasMass () const  { return problem_.hasMass(); }

    //! evaluate mass matrix factor 
    void mass(const Entity& en,
              const VolumeQuadratureType& quad,
              const int quadPoint,
              MassFactorType& m)
    {
      assert( enQuadId_ == quad.id() );

      // call problem implementation 
      problem_.mass(en, time_, quad.point(quadPoint),
                    valuesVec_[ quadPoint ],
                    m);
    }
    
  private:
    DGDiscreteModelCaller(const DGDiscreteModelCaller&);
    DGDiscreteModelCaller& operator=(const DGDiscreteModelCaller&);

  protected:
    DiscreteModelType& problem_;

  protected:  
    using BaseType::data_;
    using BaseType::time_;

    std::vector< RangeTupleType > valuesEnVec_;
    std::vector< RangeTupleType > valuesNbVec_;

    std::vector< RangeTupleType > valuesVec_;
    std::vector< JacobianRangeTupleType > jacobianVec_;

#ifndef NDEBUG
    size_t enQuadId_; 
    size_t faceQuadId_; 
#endif
  };

}

#endif
