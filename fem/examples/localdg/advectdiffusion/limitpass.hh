#ifndef DUNE_LIMITPASS_HH
#define DUNE_LIMITPASS_HH

#include <pass/pass.hh>
#include <pass/selection.hh>
#include <pass/discretemodel.hh>
#include <pass/modelcaller.hh>

#include <misc/timeutility.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

#include <vector>

namespace Dune {
  //! Concrete implementation of Pass for Limiting.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LimitDGPass :
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
    enum { dimRange = DiscreteFunctionSpaceType::DimRange,
	   dimDomain = DiscreteFunctionSpaceType::DimDomain};
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LimitDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                DiscreteFunctionSpaceType& spc) :
      BaseType(pass, spc),
      caller_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      fMat_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0),
      source_(0.0),
      identity_(1.0),
      grads_(0.0),
      diffVar_(),
      twistUtil_(spc.grid())
    {}
    
    //! Destructor
    virtual ~LimitDGPass() {
      //delete caller_;
    }
    
  private:
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      
      dest_->clear();
      
      caller_.setArgument(*arg_);
      
    }
    
    //! Some management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      caller_.finalize();
      DestinationType* U=const_cast<DestinationType*>(Element<0>::get(*arg_));
      U->assign(dest);
    }
    
    //! Perform the limitation on all elements.
    virtual void applyLocal(EntityType& en) const
    {
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType 
	BaseFunctionSetType;
      
      //- statements
      caller_.setEntity(en);
      LocalFunctionType uEn = Element<0>::get(*arg_)->localFunction(en);
      LocalFunctionType updEn = dest_->localFunction(en);
      GeometryType geom = en.geometry().type();
      const IndexSetType& iset = spc_.indexSet();
      const BaseFunctionSetType& bsetEn = spc_.getBaseFunctionSet(en);
      DofConversionUtility< PointBased > dofConversion(dimRange);
      int numBasis=updEn.numDofs()/dimRange;
      // ************************************************
      RangeType lambda(1.),ubar;
      for (int r=0;r<dimRange;r++) {
	int dofIdx = dofConversion.combinedDof(0,r);
	updEn[dofIdx]=uEn[dofIdx];
	ubar[r]=updEn[dofIdx];
      }
      if (ubar[0]<1e-4) {
	for (int r=0;r<dimRange;r++) {
	  for (int i=0;i<numBasis;i++) {
	    int dofIdx = dofConversion.combinedDof(i,r);
	    updEn[dofIdx] = 0.;
	  }
	}
      } else {
	bool modified=false;
	RangeType U;
	DomainType corner[3];
	corner[0][0] = 0.;
	corner[0][1] = 0.;
	corner[1][0] = 1.;
	corner[1][1] = 0.;
	corner[2][0] = 0.;
	corner[2][1] = 1.;
	for (int l=0;l<3;l++) {
	  uEn.evaluateLocal(en, corner[l], U);
 	  if (U[0]<0.) {
	    double lambdal=-ubar[0]/(U[0]-ubar[0]);
	    if (lambdal<0. || lambdal>1.) {
	      std::cerr << "LAMBDA: El(" << l << ") " 
			<< lambdal << " " << ubar[0] << " " 
			<< U[0] << " " << U[0]-ubar[0] << std::endl;
	    }
	    lambda[0]=(lambda[0]<lambdal)?lambda[0]:lambdal;
	    modified=true;
	  }
	}
	for (int i=1;i<numBasis;i++) {
	  int dofIdx0 = dofConversion.combinedDof(i,0);
	  updEn[dofIdx0] = lambda[0]*uEn[dofIdx0];
	  for (int r=1;r<dimRange;r++) {
	    int dofIdx = dofConversion.combinedDof(i,r);
	    if (!modified) 
	      updEn[dofIdx] = uEn[dofIdx];
	    else 
	      updEn[dofIdx] = updEn[dofIdx0]*ubar[r]/ubar[0];     
	  }
	}
      }
      for (int r=0;r<dimRange;r++) {
	for (int i=0;i<numBasis;i++) {
	  int dofIdx = dofConversion.combinedDof(i,r);
	  assert(uEn[dofIdx]==uEn[dofIdx]);
	  assert(updEn[dofIdx]==updEn[dofIdx]);
	}
      }
    }
    
  private:
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);
    
  private:

    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType source_;
    mutable RangeType identity_;
    mutable DomainType grads_;
    FieldVector<int, 0> diffVar_;
    TwistUtility<GridType> twistUtil_;
  };

  }
#endif
