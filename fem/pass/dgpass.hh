#ifndef DUNE_DGPASS_HH
#define DUNE_DGPASS_HH

#include "pass.hh"
#include "selection.hh"
#include "problem.hh"
#include "problemcaller.hh"

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  template <class ProblemImp, class DGPassTraits, class PreviousPassImp>
  class LocalDGPass :
    public LocalPass<DGPassTraits, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    // Base class
    typedef LocalPass<DGPassTraits, PreviousPassImp> BaseType;

    // Repetition of template arguments
    typedef ProblemImp ProblemType;
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename BaseType::SpaceType SpaceType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types extracted from the discrete function space type
    typedef typename SpaceType::GridType GridType;
    typedef typename SpaceType::DomainType DomainType;
    typedef typename SpaceType::RangeType RangeType;
    typedef typename SpaceType::JacobianRangeType JacobianRangeType;

    // Types extracted from the underlying grids
    typedef typename GridType::Traits::IntersectionIterator IntersectionIterator;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Types from the traits
    typedef typename DGPassTraits::DestinationType DestinationType;
    typedef typename DGPassTraits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DGPassTraits::FaceQuadratureType FaceQuadratureType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename ProblemType::SelectorType SelectorType;
    typedef ProblemCaller<
      ProblemType, ArgumentType, SelectorType> ProblemCallerType;

    // Range of the destination
    enum { dimRange = SpaceType::DimRange };
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LocalDGPass(ProblemType& problem, PreviousPassType& pass, SpaceType& spc) :
      BaseType(pass, spc),
      problem_(problem),
      caller_(0),
      arg_(0),
      dest_(0),
      spc_(spc),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      tmp_(0.0),
      grads_(0.0)
    {}
   

    virtual ~LocalDGPass() {
      //delete caller_;
    }
  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;

      if (caller_) {
        caller_->setArgument(*arg_);
      }
      else {
        // * Move this initialisation garbage into the ProblemCaller
        caller_ = new ProblemCallerType(problem_, *arg_);
      }
    }

    //! Nothing to do here, really.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {}

    //! Perform the (volume and surface) integration on all elements.
    virtual void applyLocal(EntityType& en) const
    {
      std::cout << "DGLocalPass::applyLocal()" << std::endl;
      // * A first try for an implementation - far from working...
      /*
      //- Initialise quadratures and stuff
      caller_->setEntity(en);
      LocalFunctionType updEn = dest.localFunction(en);
      GeometryType geom = en.geometry().type();
      
      VolumeQuadratureType volQuad(geom);

      double vol = volumeElement(en, volQuad);
      double vol_1 = 1.0/vol;
      double dtLocal;

      // Volumetric integral part
      for (int l = 0; l < volQuad.nop(); ++l) {
        caller_->analyticalFlux(en, volQuad, l, "result");
        
        // * add evaluation of source term here as well
        // caller_->source(...);

        for (int k = 0; k < spc_.numberOfBaseFunctions(); ++k) {
          // * to be replaced by something more standard
          const std::vector<BaseJRange>& gradients =
            spc_.getBaseFunctionSet(geom).gradients(k, volQuad);
          
          tmp_ = 0.0;
          grads_ = 0.0;
          // * umv, umtv right here?
          en.geometry().jacobianInverseTransposed
            (volQuad.point(l)).umv(gradients[l][0], grads_);
          fMat_.umv(grads_, tmp_);

          // Update vector
          for (int i = 0; i < dimRange; ++i) {
            // Watch out for the minus sign here!
            double update =
              tmp_[i]*volQuad.weight(l)*
              en.geometry().integrationElement(volQuad.point(l));

            updEn[i*dimRange + k] -= update;
          } // end for i (dimRange)
        } // end for k (baseFunctions)
      } // end for l (quadraturePoints)

      // Surface integral part
      IntersectionIterator endnit = en.iend();
      IntersectionIterator nit = en.ibegin();
      FaceQuadratureType faceQuad(nit.intersectionGlobal().type());
      
      for (; nit != endnit; ++nit) {
        if (nit.neighbor()) {
          if (nit.outside()->globalIndex() > en.globalIndex()
              || nit.outside()->partitionType() == GhostEntity) {
            
            caller_->setNeighbor(*nit);
            LocalFunctionType updNeigh = dest.localFunction(*(nit.outside()));

            for (int l = 0; l < faceQuad.nop(); ++l) {
              // * might be improved by using quadrature points directly
              // * (how to deal with quadrature points on neighbor?)
              DomainType xGlobal =
                nit.intersectionGlobal().global(faceQuad.point(l));
              DomainType xLocalEn = 
                nit.intersectionSelfLocal().global(faceQuad.point(l));
              DomainType xLocalNeigh = 
                nit.intersectionNeighborLocal().global(faceQuad.point(l));
              
              // Evaluate flux
              // * add estimate for dtLocal again
              caller_->numericalFlux(nit, faceQuad, l, valEn_, valNeigh_);

              dtLocal = 0.1;
              dtLocal =  (dtLocal < std::numeric_limits<double>::min()) ?
                dtMin_ : vol/(dtLocal*h);
              if(dtLocal < dtMin_) dtMin_ = dtLocal;
              
              // * Assumption: all elements have same number of base functions
              for (int k = 0; k < spc_.numBaseFunctions(); ++k) {
                BaseRange baseEn;
                spc_.getBaseFunctionSet(geom).evaluate(k,
                                                       diffVar_,
                                                       xLocalEn,
                                                       baseEn);

                BaseRange baseNeigh;
                spc_.getBaseFunctionSet(geom).evaluate(k,
                                                       diffVar_,
                                                       xLocalNeigh,
                                                       baseNeigh);

                for (int i = 0; i < dimRange; ++i) {
                  //assert(dimRange == 1); // * Temporary
                  updEn[i*dimRange + k] +=
                    tmp_[i]*baseEn[0]*faceQuad.weight(l)*
                    nit.intersectionGlobal().integrationElement(faceQuad.point(l));
                  updNeigh[i*dimRange + k] -= 
                    tmp_[i]*baseNeigh[0]*faceQuad.weight(l)*
                    nit.intersectionGlobal().integrationElement(faceQuad.point(l));
                } // end loop dimRange
              } // end loop base functions
            } // end loop quadrature points
            
          } // end if outside etc.
        } // end if neighbor
        if (nit.boundary()) {
          std::cout << "More to come" << std::endl;
          assert(false);
        } // end if boundary
      }
      */
    }
    
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
    ProblemType& problem_;
    mutable ProblemCallerType* caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    SpaceType& spc_;

    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType tmp_;
    mutable DomainType grads_;
    FieldVector<int, 0> diffVar_;
  };
  
} // end namespace Dune

#endif
