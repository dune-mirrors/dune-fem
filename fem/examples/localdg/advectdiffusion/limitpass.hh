#ifndef DUNE_LIMITPASS_HH
#define DUNE_LIMITPASS_HH

#include <pass/pass.hh>
#include <pass/selection.hh>
#include <pass/discretemodel.hh>
#include <pass/modelcaller.hh>

#include <misc/timeutility.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

#include <dune/quadrature/barycenter.hh>

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
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType 
    VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType 
    FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType 
    JacobianRangeType;
    
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
      diffVar_()
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
      GeometryType geom = en.geometry().type();
      //- statements
      caller_.setEntity(en);
      // get U on entity
      LocalFunctionType uEn = Element<0>::get(*arg_)->localFunction(en);
      // get local funnction for limited values
      LocalFunctionType limitEn = dest_->localFunction(en);
      BaryCenterQuad<double,DomainType,0> center0(geom);
      BaryCenterQuad<double,DomainType,1> center1(geom);
      const IndexSetType& iset = spc_.indexSet();
      const BaseFunctionSetType& bsetEn = limitEn.getBaseFunctionSet();
      Dune::DofConversionUtility< PointBased > dofConversion(dimRange);

      // number of sclar base functions
      int numBasis=limitEn.numDofs()/dimRange;
      
      double areae = volumeElement(en);

      // get value of U in barycenter
      RangeType centerUe(0);
      uEn.evaluateLocal(en,center0.point(0),centerUe);

      FieldVector<bool,dimRange> limit(false);

      IntersectionIterator endnit = en.iend();
      IntersectionIterator nit = en.ibegin();
      for (; nit != endnit; ++nit) {
	if (nit.neighbor()) {
          // get neighbor entity
	  EntityType& nb = *nit.outside(); 
	  double arean = volumeElement(nb);
	  // get local function for U on neighbor
	  LocalFunctionType uNeigh =Element<0>::get(*arg_)->localFunction(nb);
	  // get U in barycenter of neighbor
	  RangeType centerUn(0);
	  uNeigh.evaluateLocal(nb,center0.point(0),centerUn);
	  // get U on interface
	  RangeType faceUe(0);
	  RangeType faceUn(0);
	  uEn.evaluateLocal(en
			    ,*nit.intersectionSelfLocal().
			    global(center1.point(0))
			    ,faceUe);
	  uNeigh.evaluateLocal(nb
			       ,*nit.intersectionNeighborLocal().
			       global(center1.point(0))
			       ,faceUn);
	  RangeType jump = faceUn;
	  jump -= faceUe;
	  jump *= 1.0;
	  for (int r=0;r<dimRange;r++) {
	    double jumpr = jump[r];
	    if (1) 
	      limit[r] = true;
	  }
	}
      }
      for (int r=0;r<dimRange;r++) {
	if (limit[r]) {
	  int dofIdx = dofConversion.combinedDof(0,r);
	  limitEn[dofIdx] = uEn[dofIdx];
	  for (int i=1;i<numBasis;i++) {
	    int dofIdx = dofConversion.combinedDof(i,r);
	    limitEn[dofIdx] = 0.;
	  }
	}
	else {
	  for (int i=0;i<numBasis;i++) {
	    int dofIdx = dofConversion.combinedDof(i,r);
	    limitEn[dofIdx] = uEn[dofIdx];
	  }
	}
      }
    }

    
  private:
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);
    
  private:
    double volumeElement(const EntityType& en) const {
      BaryCenterQuad<double,DomainType,0> center(en);
      double result = 0;
      for (int qp = 0; qp < center.nop(); ++qp) {
        result +=
          center.weight(qp) * en.geometry().integrationElement(center.point(qp));
      }
      return result;
    }

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
  };

  }
#endif
