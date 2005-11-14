#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

#include "pass.hh"
#include "selection.hh"
#include "discretemodel.hh"
#include "modelcaller.hh"

// * needs to move
// #include "../misc/timenew.hh"
#include "../misc/timeutility.hh"

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    // Types extracted from the underlying grids
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LocalDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                DiscreteFunctionSpaceType& spc) :
      BaseType(pass, spc),
      problem_(problem),
      caller_(0),
      arg_(0),
      dest_(0),
      spc_(spc),
      dtMin_(std::numeric_limits<double>::max()),
      dtOld_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      grads_(0.0),
      time_(0),
      diffVar_()
    {}
   
    //! Destructor
    virtual ~LocalDGPass() {
      //delete caller_;
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void processTimeProvider(TimeProvider* time) {
      time_ = time;
    }

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return dtMin_;
    }

  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      dest_->clear();

      if (caller_) {
        caller_->setArgument(*arg_);
      }
      else {
        // * Move this initialisation garbage into the DiscreteModelCaller
        caller_ = new DiscreteModelCallerType(problem_, *arg_);
      }

      // time initialisation
      dtMin_ = std::numeric_limits<double>::max();
      if (time_) {
        caller_->setTime(time_->time());
      }
      else {
        caller_->setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      dtOld_ = dtMin_;
      if (time_) {
        time_->provideTimeStepEstimate(dtMin_);
      }
    }

    //! Perform the (volume and surface) integration on all elements.
    virtual void applyLocal(EntityType& en) const
    {
      //std::cout << "DGLocalPass::applyLocal()" << std::endl;
    
      //- Initialise quadratures and stuff
      caller_->setEntity(en);
      LocalFunctionType updEn = dest_->localFunction(en);
      GeometryType geom = en.geometry().type();
      
      VolumeQuadratureType volQuad(geom);

      double vol = volumeElement(en, volQuad);
      //std::cout << "Vol = " << vol << std::endl;
      double dtLocal;

      const typename DiscreteFunctionSpaceType::IndexSetType& iset = spc_.indexSet();

      // Volumetric integral part
      for (int l = 0; l < volQuad.nop(); ++l) {
        caller_->analyticalFlux(en, volQuad, l, fMat_);
        caller_->source(en, volQuad, l, source_);

        for (int k = 0; k < spc_.getBaseFunctionSet(en).numBaseFunctions(); 
             ++k) {
          JacobianRangeType gradients(0.0);
          RangeType values(0.0);
          spc_.getBaseFunctionSet(en).jacobian(k, volQuad, l, gradients);
          spc_.getBaseFunctionSet(en).eval(k, volQuad, l, values);

          // Update dof k
          // Scalar product between source contribution and base function k
          double update = source_*values;
          // Summing over all dimensions for analytical flux contribution
          for (int i = 0; i < dimRange; ++i) {
            grads_ = 0.0;
            // Scaling
            en.geometry().jacobianInverseTransposed(volQuad.point(l)).
              umv(gradients[i], grads_);
       
            // Scalar product for component i
            update += fMat_[i]*grads_;
          } // end for i (dimRange)

          // Scaling with quadrature weight and integration element
          update *= volQuad.weight(l)*
            en.geometry().integrationElement(volQuad.point(l));
          update /= vol; // * right here?

          updEn[k] += update;
        } // end for k (baseFunctions)
      } // end for l (quadraturePoints)

      // Surface integral part
      IntersectionIterator endnit = en.iend();
      IntersectionIterator nit = en.ibegin();
      FaceQuadratureType faceQuad(nit.intersectionGlobal().type());
      
      for (; nit != endnit; ++nit) {
        if (nit.neighbor()) {
          if (iset.index(*nit.outside()) > iset.index(en)
              || nit.outside()->partitionType() == GhostEntity) {
            
            caller_->setNeighbor(*nit.outside());
            LocalFunctionType updNeigh =dest_->localFunction(*(nit.outside()));

            for (int l = 0; l < faceQuad.nop(); ++l) {
              double h = 
                nit.intersectionGlobal().integrationElement(faceQuad.point(l));
              // * might be improved by using quadrature points directly
              // * (how to deal with quadrature points on neighbor?)
              DomainType xLocalEn = 
                nit.intersectionSelfLocal().global(faceQuad.point(l));
              DomainType xLocalNeigh = 
                nit.intersectionNeighborLocal().global(faceQuad.point(l));
              
              // Evaluate flux
              dtLocal = 
                caller_->numericalFlux(nit, faceQuad, l, valEn_, valNeigh_);

              dtLocal =  (dtLocal < std::numeric_limits<double>::min()) ?
                dtMin_ : vol/(dtLocal*h);
              if(dtLocal < dtMin_) dtMin_ = dtLocal;
              
              // * Assumption: all elements have same number of base functions
              for (int k = 0;
                   k < spc_.getBaseFunctionSet(en).numBaseFunctions(); ++k) {
                spc_.getBaseFunctionSet(en).evaluate(k,
                                                     diffVar_,
                                                     xLocalEn,
                                                     baseEn_);

                spc_.getBaseFunctionSet(en).evaluate(k,
                                                     diffVar_,
                                                     xLocalNeigh,
                                                     baseNeigh_);

                // * Warning: division by vol only correct for orthonormal bases!!!!!!!!!!!

                updEn[k] -=
                  (valEn_*baseEn_)*faceQuad.weight(l)/vol;
                
                updNeigh[k] += 
                  (valNeigh_*baseNeigh_)*faceQuad.weight(l)/vol;
              } // end loop base functions
            } // end loop quadrature points
            
          } // end if outside etc.
        } // end if neighbor

        if (nit.boundary()) {
          for (int l = 0; l < faceQuad.nop(); ++l) {
            double integrationElement =
              nit.intersectionGlobal().integrationElement(faceQuad.point(l));
              
            DomainType xLocalEn = 
              nit.intersectionSelfLocal().global(faceQuad.point(l));
            DomainType xGlobal =
              nit.intersectionGlobal().global(faceQuad.point(l));
            
            dtLocal = 
              caller_->boundaryFlux(nit, faceQuad, l, source_);
            dtLocal = (dtLocal < std::numeric_limits<double>::min()) ?
              dtMin_ : vol/(dtLocal*integrationElement);
            if (dtLocal < dtMin_) dtMin_ = dtLocal;
                    
            for (int k = 0; 
                 k < spc_.getBaseFunctionSet(en).numBaseFunctions(); ++k) {
              spc_.getBaseFunctionSet(en).evaluate(k,
                                                   diffVar_,
                                                   xLocalEn,
                                                   baseEn_);
              updEn[k] -= (source_*baseEn_)*faceQuad.weight(l)/vol;
            }
          }
        } // end if boundary
      }
     
    }
    
  private:
    LocalDGPass();
    LocalDGPass(const LocalDGPass&);
    LocalDGPass& operator=(const LocalDGPass&);

  private:
    double volumeElement(const EntityType& en,
                         const VolumeQuadratureType& quad) const {
      double result = 0.0;
      for (int qp = 0; qp < quad.nop(); ++qp) {
        result += 
          quad.weight(qp) * en.geometry().integrationElement(quad.point(qp));
      }
      return result;
    }
    
  private:
    DiscreteModelType& problem_;
    mutable DiscreteModelCallerType* caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
    mutable double dtMin_;
    mutable double dtOld_;

    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable DomainType grads_;
    TimeProvider* time_;
    FieldVector<int, 0> diffVar_;
  };
  
} // end namespace Dune

#endif
