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
#include <dune/fem/misc/timeprovider.hh>

#include <dune/fem/pass/ldgflux.hh>

#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>

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
    public DiscreteModelDefaultWithInsideOutSide< LimiterDefaultTraits<GlobalTraitsImp,Model> > 
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
      : model_(mod) , velocity_(0) {}

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
        return it.integrationOuterNormal( x ).two_norm();
      }
      else 
      {
        gLeft = gRight = 0.0;
        return 0.0;
      }
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
      model_.velocity(this->inside(),time,it.intersectionSelfLocal().global(x),
                      argULeft,velocity_);
      return ((it.outerNormal(x) * velocity_) < 0);
    }

  protected:
    const Model& model_;
    mutable DomainType velocity_;
  };


  /** \brief Concrete implementation of Pass for Limiting.
      The implemented Shock detection is described in detail in: 
        L. Krivodonova and J. Xin and J.-F. Remacle and N. Chevaugeon and J. E. Flaherty
        Shock detection and limiting with discontinuous Galerkin methods for hyperbolic conservation laws.
        Appl. Numer. Math., 48(3-4), pages 323-338, 2004.
    
      Link to paper:
        http://www.scorec.rpi.edu/REPORTS/2003-3.pdf

      Limiting is done by simply setting the polynomial order to zero.
  */
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

    typedef LeafGridPart< GridType > ContGridPartType; 

    typedef LagrangeDiscreteFunctionSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , ContGridPartType, 1 > LagrangeSpaceType;

    typedef DiscontinuousGalerkinSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , GridPartType, 0 > DG0SpaceType;

    typedef AdaptiveDiscreteFunction< LagrangeSpaceType > P1FunctionType; 
    typedef AdaptiveDiscreteFunction< DG0SpaceType > DG0FunctionType; 
    
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
      //contGridPart_( const_cast<GridType&> (gridPart_.grid()) ),
      //lagrangeSpace_( contGridPart_ ),
      //dg0Space_( const_cast<GridPartType&> (gridPart_)),
      //p1Function_( "p1-func" , lagrangeSpace_ ),
      //dg0Function_( "dg-0-func", dg0Space_ ),
      orderPower_( -((spc_.order()+1.0)/2.0)),
      dofConversion_(dimRange),
      faceQuadOrd_(spc_.order()),
      jump_(0),
      jump2_(0),
      gradientFlux_(1.0,10.0)
    {
      //dg0Function_.clear();
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

      /*
      // get function to limit 
      const DestinationType& U = *(Element<0>::get(*arg_));

      WeightDefault<GridPartType> weight; 
      //VtxProjectionImpl::project( U , p1Function_ , weight );
      VtxProjectionImpl::project( dg0Function_ , p1Function_ , weight );
      //p1Function_.print(std::cout);

      dg0Function_.clear();
      L2ProjectionImpl::project( U , dg0Function_ );

      {
        GrapeDataDisplay< GridType > grape(U.space().gridPart());
        grape.addData(dg0Function_);
        grape.addData(p1Function_);
        grape.dataDisplay(U);
      }
      */
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
      enum { dim = EntityType :: dimension };
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
      
      double radius = 0.0;
      const GeometryType geomType = geo.type();
      if ( geomType.isSimplex() )
      {
        const int numcorners = geo.corners();
        DomainType diff;
        for(int n=0; n<numcorners; ++n)
        {
          diff = geo[(n+1)%numcorners] - geo[n];
          radius = std::max(diff.two_norm(),radius);
        }
      }
      else if ( geomType.isCube() )
      {
        // map of used geometry points to determ side length (see
        // reference elements)
        const double hypo = SQR((geo[1][0] - geo[0][0])*0.5) + SQR((geo[2][1] - geo[0][1])*0.5);

        if( dim == 3 )
        {
          radius = sqrt(hypo + SQR((geo[4][2] - geo[0][2])* 0.5));
        } 
        else if ( dim == 2 )
        {
          // radium of circum 
          radius = sqrt(hypo);
        }
        else 
        {
          radius = std::abs((geo[1][0] - geo[0][0])*0.5);
        }
      }
      else 
      {
        DUNE_THROW(NotImplemented,"Unsupported geometry type!");
      }
            
      const double hPowPolOrder = (1.0/(geo.volume())) * pow(radius, orderPower_);// -((order+1.0)/2.0) );
      //const double hPowPolOrder = pow(4.0*M_SQRT2*radius, orderPower_);// -((order+1.0)/2.0) );
      // get value of U in barycenter

      FieldVector<bool,dimRange> limit(false);

      RangeType totaljump(0);

      // calculate circume during neighbor check 
      double circume = 0.0;

      JacobianRangeType enGrad;
      JacobianRangeType nbGrad;
      
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

            applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump,circume);
          }
          else  
          { // non-conforming case 
            typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
            NonConformingQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
            NonConformingQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

            applyLocalNeighbor(nit,faceQuadInner,faceQuadOuter,totaljump,circume);
          }
        }

        // check all neighbors 
        if (nit.boundary()) 
        {
          FaceQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
          applyBoundary(nit,faceQuadInner,totaljump,circume);
        }
      } // end intersection iterator 
       
      circume = (circume > 0.0) ? 1.0/circume : 0.0;

      // determ whether limitation is necessary  
      bool limiter = false;
      for (int r=0; r<dimRange; ++r) 
      {
        double jumpr = std::abs(totaljump[r]);
        if (jumpr*hPowPolOrder*circume > 1) 
        {
          limit[r] = true;
          limiter = true;
        }
        else
          limit[r] = false;
      }
       
      /*
      if( limiter )
      {
        typedef typename P1FunctionType :: LocalFunctionType P1LocalFunctionType;

        // get quadrature 
        VolumeQuadratureType quad(en, 2 * U.space().order() );
        const int quadNop = quad.nop();

        P1LocalFunctionType lf = p1Function_.localFunction( en ); 

        RangeType ret, tmp;
        // L2 Projection 
        {
          typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
          //! Note: BaseFunctions must be ortho-normal!!!!
          const BaseFunctionSetType& baseset = limitEn.baseFunctionSet();

          const int numDofs = limitEn.numDofs();
          for(int qP = 0; qP < quadNop ; ++qP)
          {
            lf.evaluate(quad,qP, ret);
            for(int i=0; i<numDofs; ++i)
            {
              baseset.evaluate(i,quad,qP, tmp);
              limitEn[i] += quad.weight(qP) * (ret * tmp) ;
            }
          }
        }
      }
      else 
      {
        // copy function 
        const int numDofs = limitEn.numDofs();
        for (int i=0; i<numDofs; ++i) 
        {
          limitEn[i] = uEn[i];
        }
      }
      */

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
          for (int i=1; i<numBasis; ++i) 
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

      /*
      if( limiter )
      {
        bool zeroIt = false;
        IntersectionIteratorType endnit = gridPart_.iend(en); 
        for (IntersectionIteratorType nit = gridPart_.ibegin(en); 
             nit != endnit; ++nit) 
        {
          // check all neighbors 
          if (nit.neighbor()) 
          {
            typedef TwistUtility<GridType> TwistUtilityType;
            // conforming case 
            if( TwistUtilityType::conforming(gridPart_.grid(),nit) )
            {
              FaceQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
              FaceQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

              zeroIt = limitFunc(nit,faceQuadInner,faceQuadOuter,uEn);
            }
            else  
            { // non-conforming case 
              typedef typename FaceQuadratureType :: NonConformingQuadratureType NonConformingQuadratureType;
              NonConformingQuadratureType faceQuadInner(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::INSIDE);
              NonConformingQuadratureType faceQuadOuter(gridPart_,nit,faceQuadOrd_,FaceQuadratureType::OUTSIDE);

              zeroIt = limitFunc(nit,faceQuadInner,faceQuadOuter,uEn);
            }
          }
          if ( zeroIt ) break;
        } // end intersection iterator 

        if( zeroIt )
        {
          // apply limitation 
          for (int r=0; r<dimRange; ++r) 
          {
            // dof 0 
            {
              const int dofIdx = dofConversion_.combinedDof(0,r);
              limitEn[dofIdx] = uEn[dofIdx];
            }
            // all other dofs are set to zero
            for (int i=1; i<numBasis; ++i) 
            {
              const int dofIdx = dofConversion_.combinedDof(i,r);
              limitEn[dofIdx] = 0.;
            }
          }
        }
        else 
        {
          typedef typename P1FunctionType :: LocalFunctionType P1LocalFunctionType;
            const int numDofs = limitEn.numDofs();
              for(int i=0; i<numDofs; ++i)
              {
                limitEn[i] = 0.0;
              }

          // get quadrature 
          VolumeQuadratureType quad(en, 2 * U.space().order() );
          const int quadNop = quad.nop();

          P1LocalFunctionType lf = p1Function_.localFunction( en ); 

          RangeType ret, tmp;
          // L2 Projection 
          {
            typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
            //! Note: BaseFunctions must be ortho-normal!!!!
            const BaseFunctionSetType& baseset = limitEn.baseFunctionSet();

            const int numDofs = limitEn.numDofs();
            for(int qP = 0; qP < quadNop ; ++qP)
            {
              lf.evaluate(quad,qP, ret);
              for(int i=0; i<numDofs; ++i)
              {
                baseset.evaluate(i,quad,qP, tmp);
                limitEn[i] += quad.weight(qP) * (ret * tmp) ;
              }
            }
          }
        }
      }
      else 
      {
        // copy function 
        const int numDofs = limitEn.numDofs();
        for (int i=0; i<numDofs; ++i) 
        {
          limitEn[i] = uEn[i];
        }
      }
      */
    }
    
  private:
    template <class QuadratureImp>
    void applyLocalNeighbor(IntersectionIteratorType & nit,
            const QuadratureImp & faceQuadInner,
            const QuadratureImp & faceQuadOuter,
            RangeType& totaljump,
            double& umfang) const
    {
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        double umf = caller_.numericalFlux(nit, faceQuadInner, faceQuadOuter,l,
                              jump_, jump2_);

        // if jump is not zero 
        if( umf > 0.0 )
        {
          umf *= faceQuadInner.weight(l);
          jump_ *= umf;
          umfang += umf;
          totaljump += jump_;
        }
      }
    }

    template <class QuadratureImp>
    bool limitFunc(IntersectionIteratorType & nit,
            const QuadratureImp & faceQuadInner,
            const QuadratureImp & faceQuadOuter,
            const LocalFunctionType& uEn) const 
    {
      // get function to limit 
      const DestinationType& U = *(Element<0>::get(*arg_));

      // get neighbor entity
      EntityPointerType ep = nit.outside();
      EntityType& nb = *ep; 

      // get U on entity
      const LocalFunctionType uNb = U.localFunction(nb);

      JacobianRangeType uLeft;

      //VolumeQuadratureType quad(nb,0);
      JacobianRangeType uRight;
      //uNb.jacobian( quad, 0 , uRight ); 

      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        const DomainType normal = nit.integrationOuterNormal(faceQuadInner.localPoint(l));

        uEn.jacobian( faceQuadInner, l , uLeft );
        uNb.jacobian( faceQuadOuter, l , uRight );

        double enSign = ( normal * uLeft[0] );
        double nbSign = ( normal * uRight[0] );

        if( enSign < 0.0 && nbSign > 0.0 ) return true;
        if( enSign > 0.0 && nbSign < 0.0 ) return true;
      }
      return false;
    }

    template <class QuadratureImp>
    void applyBoundary(IntersectionIteratorType & nit,
                      const QuadratureImp & faceQuadInner,
                      RangeType& totaljump,
                      double& umfang) const
    {
      const int faceQuadNop = faceQuadInner.nop();
      for(int l=0; l<faceQuadNop; ++l) 
      {
        // calculate jump 
        double umf = caller_.boundaryFlux(nit, faceQuadInner,l,jump_);

        // if jump is not zero 
        if( umf > 0.0 )
        {
          umf *= faceQuadInner.weight(l);
          jump_ *= umf;
          totaljump += jump_;
          umfang += umf; 
        }
      }
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
    //ContGridPartType contGridPart_;
    //LagrangeSpaceType lagrangeSpace_;
    //DG0SpaceType dg0Space_;
    //mutable P1FunctionType p1Function_;
    //mutable DG0FunctionType dg0Function_;
    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
    const int faceQuadOrd_;
    mutable RangeType jump_; 
    mutable RangeType jump2_;
    GradientFlux gradientFlux_;
  }; // end DGLimitPass 

} // end namespace Dune 
#endif
