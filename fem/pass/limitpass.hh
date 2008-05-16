#ifndef DUNE_LIMITERPASS_HH
#define DUNE_LIMITERPASS_HH

//- system includes 
#include <vector>

//- Dune includes 
#include <dune/common/fvector.hh>
#include <dune/fem/misc/utility.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/io/file/dgfparser/entitykey.hh>


#include <dune/fem/pass/dgpass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/selection.hh>
#include <dune/fem/solver/timeprovider.hh>

#include <dune/fem/pass/ldgflux.hh>

#include <dune/fem/space/dgspace.hh>
#include <dune/fem/space/lagrangespace.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/operator/projection/l2projection.hh>
#include <dune/fem/operator/projection/vtxprojection.hh>

//#include <dune/fem/operator/artificialdiffusion.hh>

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
    
#ifdef USE_LIMITER_AFTER
    typedef Selector<1> SelectorType;
#else
    typedef Selector<0> SelectorType;
#endif
    typedef FieldVector<double, Traits::dimDomain> DomainType;
    typedef FieldVector<double, Traits::dimDomain-1> FaceDomainType;
    typedef typename Traits::RangeType RangeType;
    typedef typename Traits::GridType GridType;
    typedef typename Traits::JacobianRangeType JacobianRangeType;
    typedef typename Traits::GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Entity EntityType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename DomainType :: field_type DomainFieldType;

    enum { dimRange = RangeType :: dimension };
    
  public:
    /** \brief default limiter discrete model */
    LimiterDefaultDiscreteModel(const Model& mod, const DomainFieldType veloEps = 1e-12) 
      : model_(mod) , velocity_(0) , veloEps_(veloEps) 
    {}

    //! \brief returns false 
    bool hasSource() const { return false; }
    //! \brief returns true 
    bool hasFlux() const   { return true;  }
    
    /** \brief numericalFlux of for limiter evaluateing the difference
         of the solution in the current integration point if we are a an
         inflow intersection.
         This is needed for the shock detection.
    */
    template <class ArgumentTuple>
    double numericalFlux(const IntersectionIteratorType& it,
                         const double time, 
                         const FaceDomainType& x,
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

    /** \brief boundaryFlux to evaluate the difference of the interior solution 
        with the boundary value. This is needed for the limiter. 
        The default returns 0, meaning that we use the interior value as ghost value. 
    */
    template <class ArgumentTuple>
    double boundaryFlux(const IntersectionIteratorType& it,
                        const double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft, 
                        RangeType& gLeft) const
    { 
      gLeft = 0 ;
      return 0.0;
    }

    void indicatorMax()
    {
    }

    /** \brief adaptation method */
    void adaptation(GridType& grid, EntityPointerType& ep, 
                    const RangeType& indicator) const 
    {
    }

  protected:
    //! returns true, if we have an inflow boundary
    template <class ArgumentTuple>
    bool checkDirection(const IntersectionIteratorType& it,
                        const double time, const FaceDomainType& x,
                        const ArgumentTuple& uLeft, 
                        const ArgumentTuple& uRight,
                        RangeType& gLeft,
                        RangeType& gRight) const 
    { 
      typedef typename ElementType<0, ArgumentTuple>::Type UType;
      const UType& argULeft = Element<0>::get(uLeft);
      model_.velocity(this->inside(),time,it.intersectionSelfLocal().global(x),
                      argULeft,velocity_);

      // calculate scalar product of normal and velocity 
      const double scalarProduct = it.outerNormal(x) * velocity_;

      // check inflow boundary 
      // in case of product zero also check 
      // (otherwise errors on problems with only diffusion)
      return (scalarProduct < 0) ? true : // inflow intersection 
               (scalarProduct > 0) ? false :  // outflow intersection
               (velocity_.two_norm() < veloEps_); // velocity zero 
    }

  protected:
    const Model& model_;
    mutable DomainType velocity_;
    const DomainFieldType veloEps_;
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
  template< class DiscreteModelImp , class PreviousPassImp , int passId  = -1 >
  class LimitDGPass :
    public LocalPass< DiscreteModelImp , PreviousPassImp , passId > 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp , PreviousPassImp , passId > BaseType;
    typedef LimitDGPass< DiscreteModelImp , PreviousPassImp , passId > ThisType;
    
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
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType 
    JacobianRangeType;
    
    // Types extracted from the underlying grids
    typedef typename GridPartType :: IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::EntityPointer EntityPointerType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;
    
    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef CombinedSelector< ThisType , SelectorType > CombinedSelectorType;
    typedef DiscreteModelCaller< DiscreteModelType 
                                 , ArgumentType 
                                 , CombinedSelectorType
                               > DiscreteModelCallerType;

    typedef DofConversionUtility< PointBased > DofConversionUtilityType;

    typedef LeafGridPart< GridType > ContGridPartType; 

    typedef LagrangeDiscreteFunctionSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , ContGridPartType, 1 > LagrangeSpaceType;

    typedef DiscontinuousGalerkinSpace < typename DiscreteFunctionSpaceType
      :: FunctionSpaceType , GridPartType, 0 > DG0SpaceType;

    typedef AdaptiveDiscreteFunction< LagrangeSpaceType > P1FunctionType; 
    typedef AdaptiveDiscreteFunction< DG0SpaceType > DG0FunctionType; 

    // Range of the destination
    enum { dimRange = DiscreteFunctionSpaceType::dimRange,
     dimDomain = DiscreteFunctionSpaceType::dimDomain};
    enum { dimension = GridType :: dimension };
    typedef typename GridType :: ctype ctype; 
    typedef FieldVector<ctype, dimDomain-1> FaceDomainType;

    //! is true if grid is structured grid 
    enum { StructuredGrid = ! Capabilities::IsUnstructured<GridType>::v };

  public:
    //- Public methods
    /** \brief constructor
     * 
     *  \param  problem    Actual problem definition (see problem.hh)
     *  \param  pass       Previous pass
     *  \param  spc        Space belonging to the discrete function local to
     *                     this pass
     *  \param  paramFile  parameter file (optional) containing whether to use
     *                     TVD switch or not:
     *                     # 0 == TVB , 1 == TVD 
     *                     TVD: 0  # defaut value 
     */
    LimitDGPass(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                const std::string paramFile = "" ) :
      BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      orderPower_( -((spc_.order()+1.0)/2.0)),
      dofConversion_(dimRange),
      faceQuadOrd_(spc_.order()),
      jump_(0),
      jump2_(0),
      gradientFlux_(1.0,10.0),
      comboSet_(),
      tvdAttribute_( getTvdParameter(paramFile) )
    {
      
      {
        typedef AllGeomTypes< typename GridPartType :: IndexSetType,
                              GridType> AllGeomTypesType;
        AllGeomTypesType allGeomTypes( gridPart_.indexSet() );
        
        // add all barycenters for codim 0
        {
          const std::vector<Dune::GeometryType>& geomTypes =
            allGeomTypes.geomTypes(0);

          for(size_t i=0; i<geomTypes.size(); ++i)
          {
            typedef typename Geometry :: ctype coordType;
            enum { dim = GridType :: dimension };
            const ReferenceElement< coordType, dim > & refElem =
                   ReferenceElements< coordType, dim >::general( geomTypes[i] );

            // store local coordinates of barycenter 
            baryCenterMap_[ geomTypes[i] ] = refElem.position(0,0);
          }
        }
        
        // add all barycenters for codim 1
        {
          const std::vector<Dune::GeometryType>& geomTypes =
            allGeomTypes.geomTypes(1);

          for(size_t i=0; i<geomTypes.size(); ++i)
          {
            typedef typename Geometry :: ctype coordType;
            enum { dim = GridType :: dimension-1 };
            const ReferenceElement< coordType, dim > & refElem =
                   ReferenceElements< coordType, dim >::general( geomTypes[i] );

            // store local coordinates of barycenter 
            faceCenterMap_[ geomTypes[i] ] = refElem.position(0,0);
          }
        }
      }

      // we need the flux here 
      assert(problem.hasFlux());
    }
    
    //! Destructor
    virtual ~LimitDGPass() {
    }
    
  protected:    
    // get tvd parameter from parameter file 
    bool getTvdParameter(const std::string& paramFile) const 
    {
      int tvd = 0;
      if( paramFile != "" ) 
      {
        readParameter( paramFile , "TVD", tvd );
      }
      return (tvd == 1) ? true : false;
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
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
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

  public:    
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      dest_->clear();
      caller_.setArgument(*arg_);

      problem_.indicatorMax();
    }
    
    //! Some management.
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      caller_.finalize();
#ifndef USE_LIMITER_AFTER
      DestinationType* U = const_cast<DestinationType*>(Element<0>::get(*arg_));
      U->assign(dest);
#endif
    }

    //! Perform the limitation on all elements.
    void applyLocal(EntityType& en) const
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
      // sanity check
      assert( numBasis * dimRange == limitEn.numDofs() );

      // if a component is true, then this component has to be limited 
      FieldVector<bool,dimRange> limit(false);
      // determ whether limitation is necessary  
      bool limiter = false;

      RangeType totaljump(0);

      // calculate circume during neighbor check 
      double circume = 0.0;

      const GeometryType geomType = geo.type();
      const DomainType enBary = geo.global( baryCenterMap_[geomType] );

      double radius = 0.0;
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
            
      // calculate h factor 
      const double hPowPolOrder = (1.0/(geo.volume())) * pow(radius, orderPower_);
      //const double hPowPolOrder = pow( sqrt(geo.volume()), orderPower_);
      //const double hPowPolOrder = pow(radius, orderPower_);

      JacobianRangeType enGrad;
      JacobianRangeType nbGrad;

      const IntersectionIteratorType endnit = gridPart_.iend(en); 
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

      /*
      // get || u ||
      RangeFieldType sum = 0;
      {
        VolumeQuadratureType quadrature( en, spc_.order() * 2 );
        const int numQuadraturePoints = quadrature.nop();
        for( int qp = 0; qp < numQuadraturePoints; ++qp )
        {
          const RangeFieldType factor
            = quadrature.weight( qp )
              * geo.integrationElement( quadrature.point( qp )) / geo.volume();

          RangeType uqp;
          uEn.evaluate( quadrature[ qp ], uqp );

          sum += factor * (uqp * uqp);
        }
        sum = std::sqrt(sum);
      }
      */

      // multiply h pol ord with circume 
      const double circFactor = (circume > 0.0) ? (hPowPolOrder /( circume )) : 0.0;

      // get grid 
      GridType& grid = const_cast<GridType&> (gridPart_.grid());
      
      IntersectionIteratorType nit = gridPart_.ibegin(en); 
      if( nit == endnit ) return ;
      EntityPointerType inside = nit.inside();

      for(int r=0; r<dimRange; ++r) 
      {
        totaljump[r] = std::abs(totaljump[r]) * circFactor;
        if( totaljump[r] > 1 )
        {
          limit[r] = true;
          limiter = true;
        }
      }
      
      // call problem adaptation for setting refinement marker 
      problem_.adaptation( grid, inside, totaljump );
      
      // prepare limitEn 
      for (int r=0; r<dimRange; ++r) 
      {
        if (limit[r]) 
        {
          // set to zero
          {
            const int dofIdx = dofConversion_.combinedDof(0,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
          
          // set to zero
          for (int i=1; i<numBasis; ++i) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = 0.;
          }
        }
        else 
        {
          // copy function since we do not limit anything
          for (int i=0; i<numBasis; ++i) 
          {
            const int dofIdx = dofConversion_.combinedDof(i,r);
            limitEn[dofIdx] = uEn[dofIdx];
          }
        }
      }

      if( limiter )
      {
        // default is zero 
        typedef FieldVector< DomainType , dimRange > DeoModType; 
        DeoModType deoMod;

        std::vector< DeoModType > deoMods;
        deoMods.reserve( 20 );
        std::vector< CheckType > comboVec;
        comboVec.reserve( 20 );

        std::vector< DomainType > barys;
        std::vector< RangeType >  nbVals;

        barys.reserve( 6 );
        nbVals.reserve( 6 );

        RangeType enVal;

        // evaluate function  
        evalAverage(en, uEn, enVal );

        // loop over all neighbors 
        for ( ; nit != endnit; ++nit) 
        {
          // check all neighbors 
          if (nit.neighbor()) 
          {
            EntityPointerType ep = nit.outside();
            EntityType& nb = *ep;
            
            // get U on entity
            const LocalFunctionType uNb = U.localFunction(nb);

            RangeType nbVal;

            // evaluate function  
            evalAverage(nb, uNb, nbVal );

            // calculate difference 
            nbVal -= enVal;

            // store value 
            nbVals.push_back(nbVal);

            const Geometry& nbGeo = nb.geometry();
            
            DomainType nbBary = nbGeo.global( baryCenterMap_[ nbGeo.type() ] );
            // calculate difference 
            nbBary -= enBary;

            // store value 
            barys.push_back(nbBary);

          } // end neighbor 

          // use ghost cell for limiting 
          if( nit.boundary() )
          {
            RangeType nbVal;
            RangeType avg(0);

            // get quadrature with lowest order  
            FaceQuadratureType faceQuadInner(gridPart_,nit, 0 ,FaceQuadratureType::INSIDE);
            
            const int quadNop = faceQuadInner.nop();
            for(int l=0; l<quadNop; ++l) 
            {
              // assume same value on ghost element
              caller_.boundaryFlux(nit, faceQuadInner, l, nbVal);
              avg += nbVal;
            }
            avg *= 1.0/(RangeFieldType) quadNop;

            // store value 
            nbVals.push_back( avg );

            typedef typename IntersectionIteratorType :: Geometry LocalGeometryType;
            const LocalGeometryType& interGeo = nit.intersectionGlobal();
            DomainType tmp (interGeo[0]);
            tmp -= enBary;
            tmp *= 2.0;

            const FaceDomainType& faceMid = faceCenterMap_[ interGeo.type() ];

            DomainType nbBary = nit.unitOuterNormal(faceMid);
            const double factor = nbBary * tmp;

            // calculate new value 
            nbBary *= factor;

            // store value 
            barys.push_back(nbBary);
          }
          
        } // end intersection iterator 

        // get number of found neighbors 
        const size_t neighbors = nbVals.size();

        // check maxima and minima 
        {
          bool allNegative = true;
          bool allPositive = true;
          // check whether all differences are positive or negative 
          for(size_t n=0; n<neighbors; ++n) 
          {
            for(int r=0; r<dimRange; ++r)
            {
              if(limit[r] ) 
              {
                if( nbVals[n][r] >= 0 ) allNegative = false;
                if( nbVals[n][r] <= 0 ) allPositive = false;
              }
            }
          }
          
          // stay with the p-zero version 
          // in case of local maximum or minimum 
          if( allNegative || allPositive ) return ;
        }

        typedef FieldMatrix<RangeFieldType, dim, dim> MatrixType;
        MatrixType matrix(0);
        MatrixType inverse(0);
        DomainType rhs(0);
        
        // create combination set 
        setupComboSet( neighbors , geomType );

        // calculate linear functions 
        // D(x) = U_i + D_i * (x - w_i)
        typedef typename ComboSetType :: iterator iterator; 
        const iterator endit = comboSet_.end();
        for(iterator it = comboSet_.begin(); it != endit; ++it) 
        {
          // get tuple of numbers 
          const KeyType& v = (*it).first; 

          // setup matrix 
          for(int i=0; i<dim; ++i) 
          {
            matrix[i] = barys[ v[i] ];
          }

          // create new instance of limiter coefficients  
          DeoModType dM;

          // invert matrix 
          const RangeFieldType det = 
            FMatrixHelp :: invertMatrix(matrix,inverse);
          
          // if matrix is regular 
          if( std::abs( det ) > 0 )
          {
            // calculate D
            for(int r=0; r<dimRange; ++r)
            {
              for(int i=0; i<dim; ++i) 
              {
                rhs[i] = nbVals[ v[i] ][r];
              }

              DomainType& D = dM[r];

              // set to zero 
              D = 0;
              // get solution 
              inverse.umv( rhs, D );
            }

            // store linear function 
            deoMods.push_back( dM );
            comboVec.push_back( (*it).second );
          }
          else 
          {
            // apply least square by adding another point 
            // this should make the linear system solvable 
             
            // creare vector with size = dim+1
            // the first dim components are equal to v 
            std::vector<int> nV( dim+1 );
            for(int i=0; i<dim; ++i) nV[i] = v[i];

            // take first point of list of pionts to check 
            CheckType check ( (*it).second );
            assert( check.size() > 0 );

            // get check iterator 
            typedef typename CheckType :: iterator CheckIteratorType; 
            const CheckIteratorType checkEnd = check.end();

            // take first entry as start value 
            CheckIteratorType checkIt = check.begin();
            
            // matrix should be regular now 
            if( checkIt == checkEnd )
            {
              // should work for 1 or 2 otherwise error 
              DUNE_THROW(InvalidStateException,"Check vector has no entries!");
            }

            // assign last element 
            nV[dim] = *checkIt ;
            
            // apply least square again 
            bool matrixSingular = 
                  applyLeastSquare( barys,
                                    nbVals,
                                    nV,
                                    dM,
                                    matrix,
                                    inverse,
                                    rhs );

            // test all other number to make matrix regular 
            while ( matrixSingular ) 
            {
              // go to next check element
              ++checkIt;
              
              // matrix should be regular now 
              if( checkIt == checkEnd )
              {
                // should work for 1 or 2 otherwise error 
                DUNE_THROW(InvalidStateException,"Matrix singular in Limiter!");
              }

              // assign last element 
              nV[dim] = *checkIt ;

              matrixSingular = 
                  applyLeastSquare( barys,
                                    nbVals,
                                    nV,
                                    dM,
                                    matrix,
                                    inverse,
                                    rhs );
            }

            // check for structured grids 
            if( StructuredGrid ) 
            {
              CheckIteratorType checkBegin = check.begin();
              // swap entries in case the found is not the first 
              if( checkIt != checkBegin )
              {
                int swap = *checkBegin; 
                *checkBegin = *checkIt; 
                *checkIt = swap;
                checkIt = check.begin();

                // apply changes to set to have 
                // only one iteration next time 
                const_cast<CheckType&> ((*it).second) = check;
              }
            }

            // erase the element that was used for 
            // applying least square from the list 
            check.erase( checkIt );

            // store linear function
            deoMods.push_back( dM );
            // store vector with points to check 
            comboVec.push_back( check );
          }
        } // end solving 

        // Limiting 
        {
          const size_t numFunctions = deoMods.size();
          // for all functions check with all values 
          for(size_t j=0; j<numFunctions; ++j) 
          {
            const std::vector<int> & v = comboVec[j];

            // loop over dimRange 
            for(int r=0; r<dimRange; ++r) 
            {
              RangeFieldType minimalFactor = 1.0;
              DomainType& D = deoMods[j][r];
              
              const size_t vSize = v.size();
              for(size_t m=0; m<vSize; ++m)
              {
                // get current number of entry 
                const size_t k = v[m];
                const DomainType& omega = barys[k];

                // evaluate values for limiter function 
                const RangeFieldType g = D * omega;
                const RangeFieldType d = nbVals[k][r];
                
                // if product smaller than zero set localFactor to zero 
                // otherwise d/g until 1 
                RangeFieldType localFactor = 
                      ( (g * d) < 0 ) ? 0 : 
                      ( (std::abs(g) > std::abs(d)) ) ? (d/g) : 1;
                
                // check rounding errors and in case do nothing  
                if( std::abs(g) * std::abs(d) < 1e-15 ) localFactor = 1;

                // take minimum 
                minimalFactor = std::min( localFactor , minimalFactor );
              }

              // scale linear function 
              D *= minimalFactor;
            }
          }
          
          // take maximum of limited functions 
          {
            RangeType max (0);
            std::vector< size_t > number(dimRange,0);
            for(int r=0; r<dimRange; ++r) 
            {
              max[r] = deoMods[0][r].two_norm();
            }

            for(size_t l=1; l<numFunctions; ++l) 
            {
              for(int r=0; r<dimRange; ++r) 
              {
                RangeFieldType D_abs = deoMods[l][r].two_norm();
                if( D_abs > max[r] ) 
                {
                  number[r] = l;
                  max[r] = D_abs;
                }
              }
            }
            
            for(int r=0; r<dimRange; ++r) 
            {
              deoMod[r] = deoMods[ number[r] ][r];
            }
          }

          // apply TVD switch to final function 
          if ( tvdAttribute_ )
          {
            tvdSwitch( geo, enBary, nbVals, deoMod );
          } 

          // L2 Projection 
          {
            // set zero dof to zero
            for(int r=0; r<dimRange; ++r) 
            {
              if( limit[r] ) 
              {
                const int dofIdx = dofConversion_.combinedDof(0,r);
                limitEn[dofIdx] = 0;
              }
            }

            RangeType ret, phi;
            // get quadrature for barycenter 
            VolumeQuadratureType quad( en, spc_.order() + 1 );

            typedef typename DiscreteFunctionSpaceType :: BaseFunctionSetType  BaseFunctionSetType;
            //! Note: BaseFunctions must be ortho-normal!!!!
            const BaseFunctionSetType& baseset = limitEn.baseFunctionSet();
            const int quadNop = quad.nop();
            for(int qP = 0; qP < quadNop ; ++qP)
            {
              // get global coordinate 
              DomainType point = geo.global( quad.point( qP ) );
              point -= enBary;

              // get quadrature weight 
              const double intel = quad.weight(qP);

              // evaluate linear function 
              for(int r=0; r<dimRange; ++r) 
              {
                ret[r] = enVal[r] + (deoMod[r] * point);
              }

              for(int r=0; r<dimRange; ++r) 
              {
                if( limit[r] )
                {
                  // project linear function 
                  for(int i=0; i<numBasis; ++i)
                  {
                    const int dofIdx = dofConversion_.combinedDof(i,r);
                    // here evaluateScalar could be used 
                    baseset.evaluate(dofIdx, quad[qP] , phi);
                    limitEn[dofIdx] += intel * (ret * phi) ;
                  }
                }
              }
            }
          }
        }
      } //end if limiter 
    }
    
  private:
    // make limiting scheme TVD 
    template <class FunctionType, class NbVectorType> 
    void tvdSwitch(const Geometry& geo, 
                   const DomainType& enBary,
                   const NbVectorType& nbVals,
                   FunctionType& deoMod) const 
    {
      const size_t neighbors = nbVals.size();
      assert( neighbors > 0 );

      RangeType uMax;
      RangeType uMin;
      const RangeType& nb0 = nbVals[0];
      // get initial maximum difference 
      for(int r=0; r<dimRange; ++r) 
      {
        uMax[r] = nb0[r];
        uMin[r] = nb0[r];
      }
      
      // get maxima and minima  
      for(size_t m=1; m<neighbors; ++m)
      {
        const RangeType& tmp = nbVals[m];
        for(int r=0; r<dimRange; ++r) 
        {
          if( tmp[r] > uMax[r] ) uMax[r] = tmp[r];
          if( tmp[r] < uMin[r] ) uMin[r] = tmp[r];
        }
      }
      
      // now take suitable scale 
      RangeType mini (1);
      const int corners = geo.corners();
      for(int i=0; i<corners; ++i) 
      {
        // get global coordinate 
        DomainType point = geo[i];
        point -= enBary;
        
        for(int r=0; r<dimRange; ++r) 
        {
          DomainType& D = deoMod[r];
          const RangeFieldType value = std::abs(D * point); 
          if( value > 0 ) 
          {
            {
              const RangeFieldType factor = std::abs(uMax[r] / value);
              mini[r] = std::min( mini[r] , factor );
            }
            {
              const RangeFieldType factor = std::abs(uMin[r] / value);
              mini[r] = std::min( mini[r] , factor );
            }
          }
        }
      }

      //std::cout << "Mini = " << mini << "\n";
      // scale linear function 
      for(int r=0; r<dimRange; ++r) 
      {
        deoMod[r] *= mini[r];
      }
    }

    // evaluate average of local function lf on entity en 
    void evalAverage(const EntityType& en, 
                     const LocalFunctionType& lf,
                     RangeType& val) const 
    {
      // get quadrature for barycenter 
      VolumeQuadratureType quad( en, spc_.order() );

      RangeType tmp; 

      // set value to zero
      val = 0;

      const int quadNop = quad.nop();
      for(int qp=0; qp<quadNop; ++qp) 
      {
        // evalaute function  
        lf.evaluate( quad[qp] , tmp );
        val += tmp;
      }

      assert( quadNop > 0 );
      // devide by number of evaluation points 
      val *= 1.0/((RangeFieldType) quadNop);
    }

    // matrix = A^T * A 
    template <class NewMatrixType, class MatrixType> 
    void multiply_AT_A(const NewMatrixType& A, MatrixType& matrix) const
    {
      assert( (int) MatrixType :: rows == (int) NewMatrixType :: cols );
      typedef typename MatrixType :: field_type value_type;

      for(int row=0; row< MatrixType :: rows; ++row)
      {
        for(int col=0; col< MatrixType :: cols; ++col)
        {
          value_type sum = 0;
          for(int k=0; k<NewMatrixType :: rows;  ++k)
          {
            sum += A[k][row] * A[k][col];
          }
          matrix[row][col] = sum;
        }
      }
    }

    // solve system by applying least square method 
    template <class BaryVectorType, class NbValVectorType,
              class FunctionType, class MatrixType, class VectorType> 
    bool applyLeastSquare(const BaryVectorType& barys,
                          const NbValVectorType& nbVals,
                          const std::vector<int> &nV,
                          FunctionType& dM, 
                          MatrixType& matrix, 
                          MatrixType& inverse, 
                          VectorType& rhs) const 
    {
      // dimension 
      enum { dim = dimension };
      // new dimension is dim + 1 
      enum { newDim = dim + 1 };

      // apply least square by adding another point 
      // this should make the linear system solvable 
       
      // need new matrix type containing one row more  
      typedef FieldMatrix<RangeFieldType, newDim , dim> NewMatrixType;
      typedef FieldVector<RangeFieldType, newDim > NewVectorType;

      // new matrix 
      NewMatrixType A ;

      assert( (int) nV.size() == newDim );

      // create matrix 
      for(int k=0; k<newDim; ++k) 
      {
        A[k] = barys[ nV[k] ];
      }

      // matrix = A^T * A 
      multiply_AT_A(A, matrix);

      // invert matrix 
      const RangeFieldType det = FMatrixHelp :: invertMatrix(matrix,inverse);

      if( std::abs( det ) > 0 )
      {
        // need new right hand side 
        NewVectorType newRhs;
        
        // calculate D
        for(int r=0; r<dimRange; ++r)
        {
          // get right hand side 
          for(int i=0; i<newDim; ++i) 
          {
            newRhs[i] = nbVals[ nV[i] ][r];
          }
          
          // convert newRhs to matrix 
          rhs = 0;
          A.umtv(newRhs, rhs);
          
          DomainType& D = dM[r];

          // set to zero 
          D = 0;
          // get solution 
          inverse.umv( rhs, D );
        }
        
        // return false because matrix is invertable 
        return false;
      }

      // matrix is still singular 
      return true;
    }

    // fill combination vector recursive 
    template <class SetType, int dim> 
    struct FillVector
    {
      static void fill(const int neighbors, 
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v) 
      {
        assert( (int) v.size() > dim );
        for(int n = start; n<neighbors; ++n)
        {
          v[dim] = n;
          FillVector<SetType, dim+1> :: fill( neighbors, 
                                              n + 1,
                                              comboSet,
                                              v);
        }
      } 
                       
    };
    
    // termination of fill combination vector 
    template <class SetType> 
    struct FillVector<SetType, dimension-1 >
    {
      enum { dim = dimension-1 };
      static void fill(const int neighbors, 
                       const int start,
                       SetType& comboSet,
                       std::vector<int>& v) 
      {
        assert( (int) v.size() == dimension );
        for(int n = start; n<neighbors; ++n)
        {
          v[dim] = n;
          comboSet.insert( VectorCompType( v, std::vector<int> () ) );
        }
      } 
    };
    
    // setup set storing combinations of linear functions 
    void setupComboSet(const int neighbors, const GeometryType& geomType) const 
    {
      // in case of non-conforming grid or grid with more than one
      // element type this has to be re-build on every element
      if( ( ! GridPartType :: conforming ) || 
          spc_.multipleGeometryTypes() || 
          (comboSet_.size() == 0) ) 
      {
        // clear set 
        comboSet_.clear();

        // maximal number of neighbors 
        std::vector<int> v(dimension,0);
        FillVector<ComboSetType, 0> :: fill(neighbors,0, comboSet_, v );

        // create set containing all numbers 
        std::set<int> constNumbers;
        typedef typename std::set<int> :: iterator el_iterator;
        for(int i = 0; i<neighbors; ++i)
        {
          constNumbers.insert(i);
        }
          
        const int checkSize = neighbors - dimension ;
        typedef typename ComboSetType :: iterator iterator; 
        const iterator endit = comboSet_.end();
        for(iterator it = comboSet_.begin(); it != endit; ++it) 
        {
          const KeyType& v = (*it).first;
          CheckType& check = const_cast<CheckType&> ((*it).second);

          // reserve memory 
          check.reserve ( checkSize );

          // get set containing all numbers 
          std::set<int> numbers (constNumbers);
          
          // remove the key values 
          for(int i=0; i<dimension; ++i) 
          {
            el_iterator el = numbers.find( v[i] );
            numbers.erase( el );
          }
          
          // generate check vector 
          el_iterator endel = numbers.end(); 
          for(el_iterator elit = numbers.begin(); 
              elit != endel; ++ elit)
          {
            check.push_back( *elit );
          }
        }
      }
    }

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

    // make private 
    LimitDGPass();
    LimitDGPass(const LimitDGPass&);
    LimitDGPass& operator=(const LimitDGPass&);
    
  private:
    mutable DiscreteModelCallerType caller_;
    const DiscreteModelType& problem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType& gridPart_;
    const double orderPower_;
    const DofConversionUtilityType dofConversion_; 
    const int faceQuadOrd_;
    mutable RangeType jump_; 
    mutable RangeType jump2_;
    GradientFlux gradientFlux_;

    mutable std::map< const GeometryType, DomainType > baryCenterMap_; 
    mutable std::map< const GeometryType, FaceDomainType > faceCenterMap_; 

    typedef DGFEntityKey<int> KeyType;
    typedef std::vector<int> CheckType;
    typedef std::pair< KeyType, CheckType > VectorCompType;
    typedef std::set< VectorCompType > ComboSetType;
    mutable ComboSetType comboSet_;

    // if true scheme is TVD 
    const bool tvdAttribute_; 
  }; // end DGLimitPass 

} // end namespace Dune 
#endif
