#ifndef DUNE_LIMITERPASS_HH
#define DUNE_LIMITERPASS_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/common/utility.hh>

#include <dune/grid/common/grid.hh>

#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/misc/timeutility.hh>

//*************************************************************
namespace Dune {  
/*! @addtogroup PassLimit
*/
  template <class GlobalTraitsImp, class Model>
  class LimiterDefaultDiscreteModel;
  
  template <class GlobalTraitsImp, class Model>
  struct LimiterDefaultTraits 
  {
    typedef typename Model::Traits ModelTraits;
    typedef typename ModelTraits::GridType GridType;

    enum { dimRange = ModelTraits::dimRange };
    enum { dimDomain = ModelTraits::dimDomain };

    typedef GlobalTraitsImp Traits; 
    typedef typename Traits::FunctionSpaceType FunctionSpaceType;

    typedef typename Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename Traits::GridPartType GridPartType;
    typedef typename Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename Traits::DestinationType DestinationType;
    typedef DestinationType DiscreteFunctionType;

    typedef typename DestinationType::DomainType DomainType;
    typedef typename DestinationType::RangeType RangeType;
    typedef typename DestinationType::JacobianRangeType JacobianRangeType;

    typedef LimiterDefaultDiscreteModel<GlobalTraitsImp,Model> DiscreteModelType;
  };
  
  
  // **********************************************
  // **********************************************
  // **********************************************
  template <class GlobalTraitsImp, class Model>
  class LimiterDefaultDiscreteModel :
    public DiscreteModelDefaultWithInsideOutSide<LimiterDefaultTraits<GlobalTraitsImp,Model> > 
  {
  public:
    typedef LimiterDefaultTraits<GlobalTraitsImp,Model> Traits;
    
    typedef Selector<0> SelectorType;
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIterator;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    
  public:
    LimiterDefaultDiscreteModel(const Model& mod) 
      : model_(mod) , velocity_(0) , normal_(0) {}

    bool hasSource() const { return false; }
    bool hasFlux() const   { return true;  }
    
    template <class ArgumentTuple>
    double numericalFlux(IntersectionIterator& it,
                         double time, const FaceDomainType& x,
                         const ArgumentTuple& uLeft, 
                         const ArgumentTuple& uRight,
                         RangeType& gLeft,
                         RangeType& gRight) const
    { 
      
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      const UType& argURight = Element<0>::get(uRight);
 
      if( checkDirection(it,time,x,uLeft,uRight,gLeft,gRight) )
      {
        gLeft  = argULeft;
        gLeft -= argURight;
        gRight = gLeft;
      }
      else 
      {
        gLeft = gRight = 0.0;
      }
      return 0.0;
    }

  protected:
    template <class ArgumentTuple>
    bool checkDirection(IntersectionIterator& it,
                        double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft, 
                        const ArgumentTuple& uRight,
                        RangeType& gLeft,
                        RangeType& gRight) const 
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      normal_ = it.outerNormal(x);
      model_.velocity(this->inside(),time,it.intersectionSelfLocal().global(x),
                      argULeft,velocity_);
      return ((normal_ * velocity_) < 0);
    }

  private:
    const Model& model_;
    mutable DomainType velocity_;
    mutable DomainType normal_;
  };


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
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType 
    JacobianRangeType;
    
    // Types extracted from the underlying grids
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;
    
    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    typedef DofConversionUtility< PointBased > DofConversionUtilityType;
    
    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::DimRange,
     dimDomain = DiscreteFunctionSpaceType::DimDomain};
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LimitDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc) :
      BaseType(pass, spc),
      caller_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      orderPower_( -((spc_.order()+1.0)/2.0)),
      dofConversion_(dimRange),
      faceQuadOrd_(spc_.order()),
      jump_(0),
      jump2_(0)
    {
      // we need the flux here 
      assert(problem.hasFlux());
    }
    
    //! Destructor
    virtual ~LimitDGPass() {
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
      DestinationType* U = const_cast<DestinationType*>(Element<0>::get(*arg_));
      U->assign(dest);
    }

    //! Perform the limitation on all elements.
    virtual void applyLocal(EntityType& en) const
    {
      // check argument is not zero
      assert( arg_ );
      
      //- statements
      // set entity to caller 
      caller_.setEntity(en);

      // get function to limit 
      const DestinationType& U = *(Element<0>::get(*arg_));

      // get U on entity
      const LocalFunctionType uEn = U.localFunction(en);

      // get local funnction for limited values
      LocalFunctionType limitEn = dest_->localFunction(en);

      // get geometry 
      const Geometry& geo = en.geometry();
      
      // number of scalar base functions
      const int numBasis = limitEn.numDofs()/dimRange;
      
      //double areae = geo.volume();
      double circume = 0.0;

      double radius = 0.0;
      
#if 0   
      const int numcorners = geo.corners();
      DomainType diff;
      for(int n=0; n<numcorners; ++n)
      {
        diff = geo[(n+1)%numcorners] - geo[n];
        double length = diff.two_norm();
        circume += length;
        radius = std::max(length,radius);

        //circume += sqrt(SQR(geo[(n+1)%numcorners][0] - geo[n][0])
        //              + SQR(geo[(n+1)%numcorners][1] - geo[n][1]));
      }
       
      // cout << "Using Alberta" << std::endl;

#else
      /*circume = 4.0*sqrt(SQR(en.geometry()[1][0] - en.geometry()[0][0]) \
  + SQR(en.geometry()[1][1] - en.geometry()[0][1]));
        
  cout << "Euclidean = " << circume << std::endl;
        
  circume = (en.geometry()[1] - en.geometry()[0]).two_norm();
        
  cout << "Two norm = " << circume << std::endl;*/
       
      // adjust for 3d 
      double ax,ay,bx,by;
        
      ax = geo[0][0];
      ay = geo[0][1];
        
      bx = geo[1][0];
      by = geo[1][1];
        
      circume = 4.0*(fabs(bx-ax)+fabs(by-ay));
        
      /*cout << "Using Yaspgrid" << std::endl;
  cout << "Point 1,x = " << en.geometry()[1][0] << std::endl;
  cout << "Point 1,y = " << en.geometry()[1][1] << std::endl;
  cout << "Point 0,x = " << en.geometry()[0][0] << std::endl;
  cout << "Point 0,y = " << en.geometry()[0][1] << std::endl;*/
        
      // cout << "Circumference = " << circume << std::endl;    
      radius = fabs(bx-ax)/M_SQRT2;   
#endif
            
      const double hPowPolOrder = pow(4.0*M_SQRT2*radius, orderPower_);// -((order+1.0)/2.0) );
      // get value of U in barycenter

      FieldVector<bool,dimRange> limit(false);

      RangeType totaljump(0);
      int counter = 0;
      
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
           nit != endnit; ++nit) 
      {
        // check all neighbors 
        if (nit.neighbor()) 
        {
          // get neighbor entity
          EntityPointerType ep = nit.outside();
          EntityType& nb = *ep; 

          // set neighbor to caller 
          caller_.setNeighbor(nb);
          
          typedef TwistUtility<GridType> TwistUtilityType;
          // conforming case 
          if( TwistUtilityType::conforming(gridPart_.grid(),nit) )
          {
            FaceQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
            FaceQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

            counter += applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump);
          }
          else  
          { // non-conforming case 
            typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
            NonConformingQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
            NonConformingQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

            counter += applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump);
          }
        }
      } // end intersection iterator 
       
      // determ whether limitation is necessary  
      for (int r=0;r<dimRange; ++r) 
      {
        double jumpr = fabs(totaljump[r])/double(counter);
        if (jumpr*hPowPolOrder > 1.) 
          limit[r] = true;
        else
          limit[r] = false;
      }
           
      // apply limitation 
      for (int r=0; r<dimRange; ++r) 
      {
        if (limit[r]) 
        {
          // dof 0 
          {
            const int dofIdx = dofConversion_.combinedDof(0,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
          // all other dofs are set to zero
          for (int i=1; i<numBasis;i++) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = 0.;
          }
        }
        else 
        {
          // copy function 
          for (int i=0; i<numBasis; ++i) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
        }
      }
    }
    
  private:
    template <class QuadratureImp>
    int applyLocalNeighbor(IntersectionIteratorType & nit,
            const QuadratureImp & faceQuadInner,
            const QuadratureImp & faceQuadOuter,
            RangeType& totaljump) const
    {
      const int faceQuadNop = faceQuadInner.nop();
      int counterInc = 0;
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,l,
                              jump_, jump2_);

        // if jump is not zero 
        if( jump_.one_norm() > 0.0 )
        {
          jump_ *= faceQuadInner.weight(l);
          totaljump += jump_;
          counterInc = 1;
        }
      }
      return counterInc; 
    }

    //! The actual computations are performed as follows. First, prepare
    //! the grid walkthrough, then call applyLocal on each entity and then
    //! call finalize.
    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // limitation only necessary if order > 0
      if( spc_.order() > 0 )
      {
        // prepare, i.e. set argument and destination 
        prepare(arg, dest);

        // dod limitation 
        IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) {
          applyLocal(*it);
        }
        // finalize
        finalize(arg, dest);
      }
      else 
      {
        // otherwise just copy 
        const DestinationType& U = *(Element<0>::get(arg));
        dest.assign(U);
      }
    }
    
    // make private 
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);
    
  private:
    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType& gridPart_;
    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
    const int faceQuadOrd_;
    mutable RangeType jump_; 
    mutable RangeType jump2_;
  }; // end DGLimitPass 

} // end namespace Dune 
#endif
