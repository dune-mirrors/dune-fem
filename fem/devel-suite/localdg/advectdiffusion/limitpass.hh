#ifndef DUNE_LIMITPASS_HH
#define DUNE_LIMITPASS_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>

//- local includes 
#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/modelcaller.hh>

#include <dune/fem/misc/timeprovider.hh>

#include <dune/fem/quadrature/elementquadrature.hh>


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
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType 
    DiscreteFunctionSpaceType;
    
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
    enum { dimRange = DiscreteFunctionSpaceType::dimRange,
     dimDomain = DiscreteFunctionSpaceType::dimDomain};
    typedef FieldVector<double, dimDomain-1> FaceDomainType;
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
      gridPart_(spc_.gridPart()),
      orderPower_( -((spc_.order()+1.0)/2.0)),
      dofConversion_(dimRange)
    {}
    
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
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
      //- statements
      caller_.setEntity(en);
      DestinationType& U = const_cast<DestinationType&> (*(Element<0>::get(*arg_)));

      // get U on entity
      LocalFunctionType uEn = U.localFunction(en);
      // get local funnction for limited values
      LocalFunctionType limitEn = dest_->localFunction(en);

      // get geometry 
      const Geometry& geo = en.geometry();
      
      // create quadrature with barycenters for 
      CachingQuadrature<GridPartType,0> center0(en,0);
      
      //const IndexSetType& iset = spc_.indexSet();
      //const BaseFunctionSetType& bsetEn = limitEn.baseFunctionSet();

      // number of scalar base functions
      const int numBasis=limitEn.numDofs()/dimRange;
      
      //double areae = geo.volume();
      double circume = 0.0;

      const int numcorners = geo.corners();
      double radius = 0.0;
      
#if 0   
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
      RangeType centerUe(0);
      uEn.evaluate(center0, 0, centerUe);

      FieldVector<bool,dimRange> limit(false);

      RangeType totaljump(0);

      int counter=0;
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
           nit != endnit; ++nit) 
      {
        if (nit.neighbor()) 
        {
          /*
          typedef TwistUtility<GridType> TwistUtilityType;
          if( TwistUtilityType::conforming(gridPart_.grid(),nit) )
          {
            typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
            FaceQuadratureType center1(gridPart_,nit,0,FaceQuadratureType::INSIDE);
          }
          */
          typedef ElementQuadrature<GridPartType,1> FaceQuadratureType;
          FaceQuadratureType center1(gridPart_,nit,0,FaceQuadratureType::INSIDE);
          FaceQuadratureType nbCenter1(gridPart_,nit,0,FaceQuadratureType::OUTSIDE);

          const DomainType normal = nit.outerNormal(center1.localPoint(0));
          //DomainType velocity(0.8);
          DomainType velocity(0.0);
          velocity[0] = 1.0;

          if (normal*velocity<0) 
          {
            // get neighbor entity
            EntityPointerType ep = nit.outside();
            const EntityType& nb = *ep; 
            // double arean = nb.geometry().volume(); 
            // get local function for U on neighbor
            
            LocalFunctionType uNeigh = U.localFunction(nb);
            // get U in barycenter of neighbor
            //RangeType centerUn(0);
            //uNeigh.evaluate(center0, 0, centerUn);

            // get U on interface
            RangeType jump(0);
            RangeType faceUn(0);

            uEn.evaluate(center1, 0, jump);
            uNeigh.evaluate(nbCenter1, 0, faceUn);
            jump -= faceUn;

            for (int r=0; r<dimRange; ++r)
              totaljump[r] += jump[r];

            ++counter;
          }
        }
      }
       
      for (int r=0;r<dimRange; ++r) 
      {
        double jumpr = fabs(totaljump[r])/double(counter);
        if (jumpr*hPowPolOrder > 1.) 
          limit[r] = true;
        else
          limit[r] = false;
      }
           
      for (int r=0;r<dimRange; ++r) 
      {
        if (limit[r]) 
        {
          int dofIdx = dofConversion_.combinedDof(0,r);
          limitEn[dofIdx] = uEn[dofIdx];
          for (int i=1;i<numBasis;i++) 
          {
            int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = 0.;
          }
        }
        else 
        {
          for (int i=0;i<numBasis; ++i) 
          {
            int dofIdx = dofConversion_.combinedDof(i,r);
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
    mutable DiscreteModelCallerType caller_;
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
    const GridPartType& gridPart_;
    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
  };

}
#endif
