#ifndef DUNE_DGPRIMALOPERATOR_HH
#define DUNE_DGPRIMALOPERATOR_HH

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/common/timer.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/fem/quadrature/caching/twistutility.hh>

//- local includes 
#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/ellipticmodelcaller.hh>

#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/space/combinedspace.hh>

#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/arrays.hh>

#include <dune/fem/misc/gridwidth.hh>

#include <dune/fem/operator/2order/dgmatrixsetup.hh>

namespace Dune {

// double feature only works in serial runs 
//#if HAVE_MPI == 0
//#define DG_DOUBLE_FEATURE 
//#endif

/*! @ingroup EllipticOperator
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/
  ////////////////////////////////////////////////////////////
  //
  //  --DGPrimalOperatorImpl 
  //
  ////////////////////////////////////////////////////////////
  /** \brief Operator assembling matrix for DG methods for elliptic problems. 
      Currently implemented are:
        - Interior Penalty
        - Baumann-Oden
        - NIPG 

        References:
          The first 4 methods can be found for example in:
            D.N. Arnold, F. Brezzi, B. Cockburn, L.D. Marini: Unified
            Analysis of Discontinuous Galerkin Methods for Elliptic
            Problems SIAM J. Num. Anal, 39  (2002), 1749-1779.
            http://www.imati.cnr.it/~marini/reports/dgm_siam.ps.gz
  */
  template <class DiscreteModelImp, 
            class PreviousPassImp, 
            class MatrixObjectTraits,
            int passId >
  class DGPrimalOperatorImpl 
    : public LocalPass<DiscreteModelImp, PreviousPassImp, passId> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp, passId> BaseType;

    typedef DGPrimalOperatorImpl<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

    //! I need to switch PreviousPassType
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType :: GridType  GridType;
    typedef typename GridType :: Traits :: LocalIdSet LocalIdSetType; 
        
    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    enum { dim = GridType :: dimension };
    
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType::polynomialOrder };

    typedef typename DiscreteModelType::Traits::ContainedSpaceType ContainedSpaceType;
    //! space of gradients of function 
    typedef CombinedSpace< ContainedSpaceType, 
                           dimGradRange,
                           PointBased > DiscreteGradientSpaceType; 

    enum { elementMassSize = DiscreteFunctionSpaceType :: localBlockSize };
    enum { massSize = DiscreteGradientSpaceType :: localBlockSize };

    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef FieldVector<double,DomainType::dimension-1> FaceDomainType;

    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename IntersectionIteratorType::Intersection IntersectionType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;

    typedef typename DiscreteModelType::SelectorType SelectorType;

    // model callers 
    typedef CombinedSelector< ThisType , SelectorType >  CombinedSelectorType;
    typedef EllipticDiscreteModelCaller< DiscreteModelType, ArgumentType,
              CombinedSelectorType> DiscreteModelCallerType;

    typedef typename GridType :: ctype ctype;
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;
    
    //my Typedefs
    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    
    typedef DiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType:: IteratorType IteratorType ;
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradJacobianRangeType;
    typedef GradJacobianRangeType GradientJacobianRangeType;
    enum { GradDimRange = GradientRangeType :: dimension };
    
    typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
     // type of temporary local function belonging to lower space 

    typedef DGMatrixTraits< MatrixObjectTraits > MyOperatorTraits;
    //! type of underlying matrix implementation 
    typedef typename MatrixObjectTraits :: template MatrixObject<
      MyOperatorTraits > :: MatrixObjectType MatrixObjectType;

    typedef typename MatrixObjectType::LocalMatrixType LocalMatrixType;
    typedef typename MatrixObjectType::MatrixType MatrixType;
    typedef typename MatrixObjectType::PreconditionMatrixType PreconditionMatrixType;
    
    typedef typename DiscreteModelType :: BoundaryIdentifierType BoundaryIdentifierType;    

    typedef typename LocalIdSetType :: IdType LocalIdType;

    typedef GradJacobianRangeType FluxRangeType; 
    typedef typename DestinationType :: LocalFunctionType SingleLFType; 

    //! singleton list , key type is const pointer to grid 
    typedef GridWidthProvider< GridType > GridWidthType;
    typedef typename GridWidthType :: ProviderType GridWidthProviderType;

    typedef typename DiscreteGradientSpaceType :: RangeType GradRangeType;

    ////////////////////////////////////////////////////
    //
    //  Coefficient and RHS caller 
    //
    ////////////////////////////////////////////////////
    template <class CallerType> 
    class CoefficientCallerTrue
    {
    public:  
      template <class QuadratureType, class CoeffType> 
      double evaluateCoefficient(CallerType& caller, 
                               EntityType& en, 
                               QuadratureType& quad, 
                               const int l,
                               CoeffType& coeff) const 
      {
        caller.evaluateCoefficient(en, quad, l, coeff );
        return coeff.infinity_norm();
      }         

      template <class CoeffType, class PsiType> 
      void applyCoefficient(const CoeffType& coeffEn, 
                            const PsiType& psi,
                            PsiType& coeffPsi) const 
      {
        for(int i=0; i<dimRange; ++i)
        {
          coeffPsi[i] = 0.0;
          coeffEn.umv(psi[i],coeffPsi[i]);
        }
      }

      template <class CoeffType, class ctype> 
      void applyCoefficient(const CoeffType& coeffEn, 
                            const FieldVector<ctype,dimDomain>& psi,
                            FieldVector<ctype,dimDomain>& coeffPsi) const 
      {
        coeffPsi = 0.0;
        coeffEn.umv(psi,coeffPsi);
      }
    };

    template <class CallerType> 
    struct CoefficientCallerFalse
    {
      template <class QuadratureType, class CoeffType> 
      double evaluateCoefficient(const CallerType& caller, 
                               const EntityType& en, 
                               const QuadratureType& quad, 
                               const int l,
                               CoeffType& coeff) const 
      {
        return 1.0;
      }         

      // just copy in default case 
      template <class CoeffType, class PsiType> 
      void applyCoefficient(const CoeffType& coeffEn, 
                            const PsiType& psi,
                            PsiType& coeffPsi) const 
      {
        coeffPsi = psi;
      }
    };


    //! if right hand side available 
    template <class CallerType> 
    class CoefficientCallerRHS 
    {
      mutable SingleLFType& singleRhs_;
      const BaseFunctionSetType& bsetEn_;
      mutable RangeType rhsval_;
      const int numDofs_ ;

    public:
      CoefficientCallerRHS(SingleLFType& singleRhs)
        : singleRhs_ ( singleRhs )
        , bsetEn_( singleRhs_.baseFunctionSet()) 
        , rhsval_ (0.0)  
        , numDofs_ ( singleRhs_.numDofs () )
      {}

      template <class QuadratureType> 
      void rightHandSide(CallerType& caller, 
                         EntityType& en, 
                         QuadratureType& volQuad, 
                         const int l, 
                         const double intel) const 
      {
        // eval rightHandSide function 
        // if empty, rhs stays 0.0
        caller.rightHandSide(en, volQuad, l, rhsval_ );

        // scale with intel 
        rhsval_ *= intel;

        //std::cout << rhsval_ << "\n";

        for (int j = 0; j < numDofs_; ++j) 
        {
          singleRhs_[j] += bsetEn_.evaluateSingle(j, volQuad[l] , rhsval_ );
        }
      }
    };

    //! no rhs 
    template <class CallerType> 
    class CoefficientCallerNoRHS 
    {
    public:  
      template <class QuadratureType> 
      void rightHandSide(const CallerType& caller, 
                         const EntityType& en, 
                         const QuadratureType& volQuad, 
                         const int l,
                         const double intel) const 
      {
      }
    };

    template <class CallerType, bool hasCoeff, bool hasRHS> 
    class CoefficientCaller : public CoefficientCallerTrue<CallerType> ,
                              public CoefficientCallerRHS<CallerType> 
    {
      typedef CoefficientCallerRHS<CallerType> BaseType;
    public:
      CoefficientCaller(SingleLFType& singleRhs) : BaseType( singleRhs ) {}
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,false,true>
      : public CoefficientCallerFalse<CallerType> ,
        public CoefficientCallerRHS<CallerType> 
    {
      typedef CoefficientCallerRHS<CallerType> BaseType;
    public:
      CoefficientCaller(SingleLFType& singleRhs) : BaseType( singleRhs ) {}
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,true,false> 
      : public CoefficientCallerTrue<CallerType>
      , public CoefficientCallerNoRHS<CallerType>
    {
      public:
    };

    template <class CallerType> 
    class CoefficientCaller<CallerType,false,false> 
      : public CoefficientCallerFalse<CallerType>
      , public CoefficientCallerNoRHS<CallerType>
    {
      public:
    };

  public:
    //- Public methods
    /**  \brief Constructor
     \param problem Actual problem definition (see problem.hh)
     \param pass Previous pass
     \param spc Space belonging to the discrete function local to this pass
     \param paramFile parameter file to read necessary parameters, if empty 
             default parameters will be applied 
    
     \note Available methods are (chosen by parameters B_{+,-}, beta
          - Interior Penalty : B_{+,-}: 0 , beta: > 0 (big)
          - Baumann-Oden     : B_{+,-}: 1 , beta: = 0       (needs polOrd > 1) 
          - NIPG             : B_{+,-}: 1 , beta: > 0      
    */         
    DGPrimalOperatorImpl(DiscreteModelType& problem, 
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                const std::string paramFile = "")
      : BaseType(pass, spc),
        caller_(problem),
        problem_(problem),
        arg_(0),
        rhs_(0),
        uh_(0),
        spc_(spc),
        gridPart_(const_cast<GridPartType &> (spc_.gridPart())),
        gradientSpace_(gridPart_),
        localIdSet_(gridPart_.grid().localIdSet()),
        gridWidth_ ( GridWidthProviderType :: getObject( &spc_.grid())),
        volumeQuadOrd_(2* spc_.order()+2),
        faceQuadOrd_(2*spc_.order()+2 ),
        matrixObj_(spc_,spc_, paramFile ),
        coeffEn_(1.0),
        coeffNb_(1.0),
        matrixAssembled_(false),
        betaFactor_(0.0),
        globalBeta_(0.0),
        beta_(0.0),
        bilinearPlus_(true),
        power_( 2*spc_.order() ),
        betaNotZero_(false),
        dtMin_ (std::numeric_limits<double>::max()),
        minLimit_(2.0*std::numeric_limits<double>::min()),
        timeDependent_( false ),
        sequence_ ( -1 )
    {
#ifndef DG_DOUBLE_FEATURE
      // we need global upwind vector to select edges 
      upwind_ = M_PI;
      if( dim > 1 ) upwind_[1] = M_LN2;
      if( dim > 2 ) upwind_[2] = M_E;
#endif
      // set to unit matrix 
      setToUnitMatrix(coeffEn_, coeffNb_);

      if( ! (spc_.order() > 0))
      {
        std::cerr << "ERROR: DG Primal operator only working for spaces with polynomial order > 0! \n";
        assert(false);
        abort();
      }
      
      bool success = true;  
      const bool output = (gridPart_.grid().comm().rank() == 0);
      // if parameter file is not empty read parameter 
      if(paramFile != "")
      {
        if( ! readParameter(paramFile,"beta",betaFactor_, output) )
          success = false;
        int bplus = 1;
        if( ! readParameter(paramFile,"B_{+,-}",bplus, output) )
          success = false;

        assert( (bplus == 0) || (bplus == 1) ); 
        bilinearPlus_ = (bplus == 0) ? false : true; 
      }

      if( ! success )
      {
        if( output ) 
        {
          std::cerr << "\nERROR: Couldn't read parameter! \n";
          std::cerr << "DGPrimalOperatorImpl -- Available Methods:\n";
          std::cerr << "Interior Penalty : B_{+,-}: 0 , beta: >  0 (big) \n";
          std::cerr << "Baumann-Oden     : B_{+,-}: 1 , beta: =  0       \n";
          std::cerr << "NIPG             : B_{+,-}: 1 , beta: >  0       \n";
        }
        exit(1);
      }

      betaNotZero_ = (std::abs(betaFactor_) > 0.0);

      if( ! betaNotZero_ )
      {
        std::cout << "DGPrimalOperatorImpl: using Baumann-Oden method!\n"; 
        if(spc_.order() < 1)
        {
          std::cerr << "WARNING: Baumann-Oden method only working for polynomial degree >= 2! \n";
          assert(false);
          exit(1);
        }
      }
      else 
      {
        if(bilinearPlus_ )
        {
          std::cout << "DGPrimalOperatorImpl: using NIPG method, beta = " << betaFactor_ << " !\n"; 
        }
        else 
        {
          std::cout << "DGPrimalOperatorImpl: using Interior Penalty method, beta = " << betaFactor_ << " !\n"; 
        }
      }
      
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );

      if(problem_.hasSource())
      {
        std::cerr << "Source for DGElliptPass not supported yet. \n";
        abort();
      }
    }

    void setToUnitMatrix(FluxRangeType& coeffEn, FluxRangeType& coeffNb) const
    {
      for( int i=0; i<FluxRangeType :: rows ; ++i) 
      {
        // set diagonal to 1 
        coeffEn[i][i] = coeffNb[i][i] = 1.0;
        // all others to zero 
        for(int j=i+1; j< FluxRangeType :: cols; ++j) 
        {
          // set default value to fMat 
          coeffEn[i][j] = coeffNb[i][j] = 0.0;
          coeffEn[j][i] = coeffNb[j][i] = 0.0;
        }
      }
    }

    //! Destructor
    virtual ~DGPrimalOperatorImpl() 
    {
    }

    void switchUpwind( double x, double y) 
    {
      upwind_[0] *= x;
      upwind_[1] *= y;
    }

    // set tau and theta 
    void enableTimeDerivative()
    {
      timeDependent_ = true ;
    }

    // set tau and theta 
    void setTauAndTheta(const double tau, const double theta)
    {
      tau_1_ = 1.0/tau;
      theta_ = theta;
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       const DestinationType &uh, 
                       DestinationType & rhs)
    {
      // store uh_ ;
      uh_ = &uh;

      computeMatrix( arg, rhs );
    }

    void buildRhs(DestinationType & rhs)
    {
      // rhs.clear();
      typedef typename DiscreteFunctionSpaceType :: IteratorType IteratorType;
      const IteratorType end = spc_.end();
      rhs_ = & rhs ;
      for(IteratorType it = spc_.begin(); it != end; ++it)
      {
        EntityType& entity = *it ; 
        if( entity.hasBoundaryIntersections() ) 
        {
          applyBoundary( *it );
        }
      }
      rhs_ = 0;
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       const bool rebuild )
    {
      if( rebuild || sequence_ != spc_.sequence() )
      {
        // time initialisation to max value 
        dtMin_ = std::numeric_limits<double>::max();

        // prepare operator  
        prepare( arg );

        // if grid has changed, then matrix structure is re-build
        matrixObj_.reserve();

        // clear matrix 
        matrixObj_.clear();
        
        // build matrix 
        const IteratorType endit = spc_.end();
        for (IteratorType it = spc_.begin(); it != endit; ++it) 
        {
          applyLocal( *it );
        }

        // build matrix and rhs 
        matrixAssembled_ = true;
        
        // create pre-condition matrix if activated 
        matrixObj_.createPreconditionMatrix();
        
        // finalize 
        finalize( );

        // store current sequence 
        sequence_ = spc_.sequence();
      }
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, 
                       DestinationType & rhs, const bool rebuild = false )
    {
      if( rebuild || sequence_ != spc_.sequence() )
      {
        //set right hand side 
        setRightHandSide( rhs );

        // compute matrix 
        computeMatrix( arg, rebuild );
      }
    }

  public:  
    //! do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      matrixObj_.multOEM( arg.leakPointer(), dest.leakPointer() ); 
    }
    
  public:
    //set right hand side 
    void setRightHandSide( DestinationType& rhs ) const 
    {
      rhs.clear ();
      rhs_ = & rhs ;
    }

    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    void prepare(const ArgumentType& arg) const
    {
      // set argument 
      arg_ = const_cast<ArgumentType*>(&arg);
      caller_.setArgument(*arg_);

      // only calculate in case of beta not zero 
      if( betaNotZero_ )
      {
        // calculate beta = O(1/h)
        globalBeta_ = betaFactor_/gridWidth_.gridWidth();
      }

      // set time to caller 
      caller_.setTime( this->time() );
    }

    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // prepare operator 
      prepare( arg );

      //set right hand side 
      setRightHandSide( dest );
    }

    //! Some timestep size management.
    void finalize() const
    {
      caller_.finalize();
      rhs_ = 0;
      uh_ = 0;
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      finalize ();
    }

    //! Estimate for the timestep size 
    double timeStepEstimateImpl() const
    {
      // factor for LDG  Discretization 
      const double p  = 2 * spc_.order() + 1;
      return dtMin_ / p;
    }

  public:
    //! return reference to function space 
    const DiscreteFunctionSpaceType & space() const 
    {
      return spc_;
    }

    //! set all data belonging to this entity to zero 
    void restrictLocal(const EntityType& father, const EntityType & son, bool firstCall ) const
    {
    }

    //! set all data belonging to this entity to zero 
    void prolongLocal(const EntityType& father, const EntityType & son, bool ) const
    {
    }

protected:
    void resizeCaches(const int numDofs) const
    {
      // resize caches 
      if(tau_.size() != numDofs) 
      { 
        tau_.resize(numDofs);
        tauNeigh_.resize(numDofs);
        phi_.resize(numDofs);
        phiNeigh_.resize(numDofs);
        psi_.resize(numDofs);
        coeffPsi_.resize(numDofs);
      }
    }

    //! assemble right hand side modifications 
    void applyBoundary(EntityType& en) const
    {
      // this method should not be called for ghost entities 
      assert( en.partitionType() != GhostEntity );
      
      // make entities known in callers
      caller_.setEntity(en);

      // create volume quadrature  
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      // get geometry
      const GeometryType & geo = en.geometry();

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );
      // resize caches 
      resizeCaches(numDofs);

      // calculate local beta factor 
      if( betaNotZero_ ) 
      {
        double betaEst = 0.0;
        if( problem_.hasCoefficient() )
        {
          CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
          betaEst = getBetaEstimate(coeffCaller, en, volQuad); 
        }
        else 
        {
          CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
          betaEst = getBetaEstimate(coeffCaller, en, volQuad); 
        }
        // calculate local beta 
        beta_ = (betaEst * globalBeta_);
      }

      assert( rhs_ );
      // local function for right hand side 
      SingleLFType singleRhs = rhs_->localFunction(en);
      
      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      const IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // neighbor volume  
        double wspeedS = 0.0;
  
        // if intersection with boundary 
        // --boundary
        if( nit.boundary() ) 
        { 
          applyLocalBoundary( *nit, 
              en, geo, volQuad, numDofs, 
              bsetEn, (LocalMatrixType* ) 0, 
              singleRhs, 
              wspeedS ); 
        } // end if boundary
      } // end intersection iterator 
    } // end applyBoundary  

    template<class QuadratureType, class CoeffCallerType> 
    double volumetricPart(EntityType& en, 
                        const GeometryType& geo,
                        QuadratureType& volQuad,
                        const CoeffCallerType& coeffCaller,
                        const BaseFunctionSetType bsetEn,
                        const int numDofs, 
                        LocalMatrixType& matrixEn) const
    {
      const int quadNop = volQuad.nop();

      double betaEst = 0.0;

      // loop over all quadrature points 
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc integration element 
        const double intel = volQuad.weight(l)
            *geo.integrationElement(volQuad.point(l));

        const JacobianInverseType& inv =
          geo.jacobianInverseTransposed(volQuad.point(l));

        ////////////////////////////////////
        // create rightHandSide
        ////////////////////////////////////
        coeffCaller.rightHandSide(caller_, en, volQuad, l, intel);
        
        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        
        // call anayltical flux of discrete model 
        betaEst = std::max(coeffCaller.evaluateCoefficient(caller_, en, volQuad, l, coeffEn_ ),
                           betaEst);

        /////////////////////////////////
        // fill element matrix 
        /////////////////////////////////
        for(int k = 0; k < numDofs; ++k)
        {
          JacobianRangeType& psi = psi_[k]; 
          JacobianRangeType& coeffPsi = coeffPsi_[k];

          // eval grad psi on reference element
          bsetEn.jacobian( k, volQuad[l], psitmp_ );
  
          // apply inverse jacobian 
          for(int i=0; i<dimRange; ++i) 
          {
            psi[i] = 0.0; 
            inv.umv(psitmp_[i], psi[i]);
          }

          // apply coefficient 
          coeffCaller.applyCoefficient(coeffEn_, psi, coeffPsi); 
        }

        // fill element matrix 
        for(int k = 0; k < numDofs; ++k)
        {
          // add diagonal entry
          {
            RangeFieldType val = 0.0;
            for (int i = 0; i <dimRange; ++i) 
            {
              val += coeffPsi_[k][i] * psi_[k][i];
            }
            val *= intel;
            matrixEn.add( k , k , val );
          }

          // add secondary diagonal
          // assume matrix is symectric
          // entry (k,j) == entry (j,k)
          for (int j = k+1; j < numDofs; ++j) 
          {
            RangeFieldType val = 0.0;
            for (int i = 0; i <dimRange; ++i) 
            {
              val += coeffPsi_[k][i] * psi_[j][i];
            }
            val *= intel;
            
            // add k,j 
            matrixEn.add( k , j , val );

            // add j,k
            matrixEn.add( j , k , val );
          }
        }
      } // end element integral 
  
      return betaEst;
    }
 
    ///////////////////////////////////////////
    //! --apply operator on entity 
    ///////////////////////////////////////////
    void applyLocal(EntityType& en) const
    {
      // this method should not be called for ghost entities 
      assert( en.partitionType() != GhostEntity );
      
      // get local element matrix 
      LocalMatrixType matrixEn = matrixObj_.localMatrix(en,en); 

      // make entities known in callers
      caller_.setEntity(en);

      // create volume quadrature  
      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      // get geometry
      const GeometryType & geo = en.geometry();

      factorFaces_ = 1;//((geo.type().isSimplex()) ? (dim+1) : 2 * dim);

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );

      double betaEst = 0.0;

      const bool rightHandSide = problem_.hasRHS() && rhs_; 

      // local function for right hand side 
      SingleLFType* singleRhsPtr = ( rhs_ ) ? 
        new SingleLFType(rhs_->localFunction(en)) :  0;
      
      // local function for right hand side 
      SingleLFType& singleRhs = *singleRhsPtr; //rhs_->localFunction(en); //rhs
      
      // resize caches 
      resizeCaches(numDofs);


      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////
      if(problem_.hasCoefficient() && rightHandSide )
      {
        CoefficientCaller<DiscreteModelCallerType,true,true> coeffCaller( singleRhs ); 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else if( problem_.hasCoefficient() )
      {
        CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else if ( rightHandSide )
      {
        CoefficientCaller<DiscreteModelCallerType,false,true> coeffCaller( singleRhs ); 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else 
      {
        CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
  
      // get beta estimate 
      beta_ = (betaEst * globalBeta_);

      const double enVolume = geo.volume();
      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      const IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // neighbor volume  
        double nbVolume = enVolume;
        double wspeedS = 0.0;

        const IntersectionType& intersection = * nit; 
  
        // if neighbor exists 
        if ( intersection.neighbor() ) 
        {
          EntityPointerType neighEp = intersection.outside();
          EntityType& nb = *neighEp;

#ifdef DG_DOUBLE_FEATURE
          // get partition type 
          const bool ghostEntity = 
            ( nb.partitionType() == GhostEntity );
          // only once per intersection or when outside is not interior 
          if( (localIdSet_.id(en) < localIdSet_.id(nb)) 
              || ghostEntity
            )
#endif
          {
            // type of TwistUtility 
            typedef TwistUtility<GridType> TwistUtilityType;
            // check conformity 
            if( nit->conforming() )
            {
              FaceQuadratureType faceQuadInner(gridPart_, intersection, faceQuadOrd_,
                                               FaceQuadratureType::INSIDE);

        
              FaceQuadratureType faceQuadOuter(gridPart_, intersection, faceQuadOrd_,
                                               FaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              nbVolume = applyLocalNeighbor( intersection, 
                              en,nb,volQuad,
                              faceQuadInner,faceQuadOuter, 
                              bsetEn,matrixEn, singleRhs, wspeedS  
#ifdef DG_DOUBLE_FEATURE
                              //, true
                            , ! ghostEntity 
#endif
                           );
            }
            else 
            {
              // we only should get here whne a non-conforming situation 
              // occurs in a non-conforming grid 
              assert( GridPartType :: conforming == false );
              
              typedef typename FaceQuadratureType :: NonConformingQuadratureType 
                NonConformingFaceQuadratureType;
              
              NonConformingFaceQuadratureType 
                nonConformingFaceQuadInner(gridPart_, intersection, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::INSIDE);
          
              NonConformingFaceQuadratureType 
                nonConformingFaceQuadOuter(gridPart_, intersection, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              nbVolume = applyLocalNeighbor(intersection,
                            en,nb,volQuad,
                            nonConformingFaceQuadInner,
                            nonConformingFaceQuadOuter, 
                            bsetEn,matrixEn, singleRhs , wspeedS 
#ifdef DG_DOUBLE_FEATURE
                            //, true
                            , ! ghostEntity 
#endif
                          );
            }
          }


        } // end if neighbor 

        // if intersection with boundary 
        // --boundary
        if( intersection.boundary() ) 
        { 
          applyLocalBoundary( intersection, 
              en, geo, volQuad, numDofs, 
              bsetEn, &matrixEn, singleRhs, 
              wspeedS ); 

        } // end if boundary

        // check timestep size 
        if ( wspeedS > minLimit_ )
        {
          const double minvolS = std::min(enVolume, nbVolume);
          const double p1 = 4.0 * (spc_.order() + 1);
          wspeedS *= p1 ;
          dtMin_ = std::min(dtMin_,minvolS/wspeedS);
        }

      } // end intersection iterator 

      // apply mass matrix to previous solution
      if( timeDependent_ )
      {
        typedef FieldMatrix<double, elementMassSize , elementMassSize > MassMatrixType; 
        MassMatrixType massEn(0);

        // mass matrix for en 
        getMassMatrix(geo,volQuad,bsetEn,numDofs,phi_,massEn);

        // apply 1/tau
        massEn *= tau_1_  ;

        //if( theta_ >=  1.0 ) 

        assert( uh_ );
        // local function for right hand side 
        const SingleLFType uhLf = uh_->localFunction(en); //rhs
      
        if( theta_ < 1.0 )
        {
          // apply matrix to old right hand side 
          multLocal( matrixEn, singleRhs, uhLf );
          matrixEn.scale( theta_ );
        }

        {
          // add to system matrix and adjust right hand side 
          for(int i=0; i<numDofs; ++i) 
          {
            for(int j=0; j<numDofs; ++j) 
            {
              singleRhs[i] += massEn[ i ][ j ] * uhLf[ j ];
              matrixEn.add( i, j, massEn[ i ][ j ] );
            }
          }
        }
      }

      // resort corresponding matrix rows for ascending numbering 
      matrixEn.resort(); 

      if( singleRhsPtr ) delete singleRhsPtr;

    } // end apply local 

    void multLocal(LocalMatrixType& matrix, SingleLFType& rhsLF, 
                   const SingleLFType& uhLf) const 
    {
      const double factor = 1.0 - theta_ ;
      typedef FieldVector<double, elementMassSize > MassVectorType; 
      MassVectorType rhs (0);

      matrix.multiplyAdd( uhLf, rhs );

      const int numDofs = rhsLF.numDofs();
      // add to system matrix and adjust right hand side 
      for(int i=0; i<numDofs; ++i) 
      {
        rhsLF[i] += factor * rhs[ i ];
      }
    }

    template <class CoeffCallerImp>
    double getBetaEstimate(CoeffCallerImp& coeffCaller, 
                           EntityType& en,
                           VolumeQuadratureType& volQuad) const 
    {
      double betaEst = 0.0;
      // loop over quadrature points 
      const int quadNop = volQuad.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // call anayltical flux of discrete model 
        betaEst = std::max(coeffCaller.evaluateCoefficient(caller_, en, volQuad, l, coeffEn_ ),
                           betaEst);
      }
      return betaEst ;
    }

    //! apply boundary integrals to matrix and right hand side 
    void applyLocalBoundary(const IntersectionType& nit, 
                            EntityType& en,
                            const GeometryType& geo, 
                            VolumeQuadratureType& volQuad,
                            const int numDofs,
                            const BaseFunctionSetType& bsetEn, 
                            LocalMatrixType* matrixEnPtr, 
                            SingleLFType& singleRhs,
                            double& wspeedS) const 
    {
      // create quadrature 
      FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                       FaceQuadratureType::INSIDE);

      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
      const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
      // get number of base functions for gradient space 

      LocalMatrixType& matrixEn = *matrixEnPtr;

      // loop over quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calculate normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol;

        const double intelFactor = faceQuadInner.weight(l);
        
        // integration element factor 
        const double intel = intelFactor * faceVol;

        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        // get boundary value 
        RangeType boundaryValue(0.0);

        // call boundary value function 
        BoundaryIdentifierType bndType = 
          caller_.boundaryValue(nit, faceQuadInner, l, boundaryValue);

        // only Dirichlet and Neumann Boundary supported right now 
        assert( bndType.isDirichletType() || bndType.isNeumannType() );

        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        assert( psi_.size() > 0 );
        JacobianRangeType& norm = psi_[0];
        if( problem_.hasCoefficient() )
        {
          // evaluate coefficient on boundary
          caller_.evaluateCoefficientBoundary(nit, faceQuadInner,l,coeffEn_);
          for(int i=0; i<dimRange; ++i)
          {
            norm[i] = 0.0;
            coeffEn_.umv(unitNormal,norm[i]);
          }
        }
        else 
        {
          for(int i=0; i<dimRange; ++i)
          {
            for(int j=0; j<dimDomain; ++j)
            {
              norm[i][j] = unitNormal[j];
            }
          }
        }

        double ldt = 0.0 ;
        // overall beta factor 
        const double facBeta = 
          factorBeta(intelFactor,faceVol, coeffEn_, coeffEn_, ldt );

        wspeedS += ldt * faceQuadInner.weight(l);

        // cache base functions evaluations
        for(int k=0; k<numDofs; ++k)
        { 
          // evaluate normal * grad phi 
          tau_[k] = bsetEn.evaluateGradientSingle(k,en, faceQuadInner[l] , norm);  
          // evaluate phi 
          bsetEn.evaluate(k,faceQuadInner[l] , phi_[k]);
        }

        // if not Babuska-Zlamal method, add boundary terms 
        {
          // only change right hand side if exists 
          if( bndType.isDirichletNonZero() && rhs_ )
          {
            // fill right hand side  
            for(int k=0; k<numDofs; ++k)
            {  
              // only valid for dim range = 1
              RangeFieldType rhsVal1 = boundaryValue[0] * tau_[k];

              rhsVal1 *= bilinIntel;
              singleRhs[k] += rhsVal1;
            }
          }

          // only on non Neumann type boundaries
          if( matrixEnPtr && bndType.isDirichletType() )
          {
            // fill matrix entries 
            for(int k=0; k<numDofs; ++k)
            {  
              for (int j = 0; j < numDofs; ++j) 
              {
                {
                  // grad w * v 
                  RangeFieldType val = tau_[j] * phi_[k][0];
                  val *= -intel;
                  matrixEn.add( k , j , val );
                }
                
                {
                  // w * grad v
                  RangeFieldType val = tau_[k] * phi_[j][0];
                  val *= bilinIntel;
                  matrixEn.add( k , j , val );
                }
              }
            }
          }
        }
            
        // dirichlet boundary values for u 
        // only change right hand side if exists 
        if(bndType.isNeumannNonZero() && rhs_ )
        {
          // fill matrix entries 
          for(int k=0; k<numDofs; ++k)
          {  
            // only valid for dim range = 1
            RangeFieldType rhsVal = boundaryValue * phi_[k];

            rhsVal *= intelFactor;
            singleRhs[k] += rhsVal;
          }
        }

        if( betaNotZero_ )
        {
          // stabilization 
          if( matrixEnPtr && bndType.isDirichletType() )
          {
            // fill matrix entries 
            for(int k=0; k<numDofs; ++k)
            {  
              for (int j = 0; j < numDofs; ++j) 
              {
                // phi_j * phi_k 
                RangeFieldType phiVal = phi_[j] * phi_[k]; 
                phiVal *= facBeta;
                matrixEn.add( k , j , phiVal );
              }
            }
            
            // dirichlet boundary values for u 
            // only change right hand side if exists 
            if(bndType.isDirichletNonZero() && rhs_ )
            {
              // fill right hand side 
              for(int k=0; k<numDofs; ++k)
              {  
                // only valid for dim range = 1
                RangeFieldType rhsVal1 = boundaryValue[0] * phi_[k];
                rhsVal1 *= facBeta;
                singleRhs[k] += rhsVal1;
              }
            } 
          }
        }
      }
    } // end applyLocalBoundary

    double compBetaK(const FluxRangeType& K) const 
    {
      double detK = K[0][0]*K[1][1]-K[1][0]*K[0][1];
      double p = (K[0][0]+K[1][1])*0.5;
      double q = p*p-detK;
      if( q < 0 && q > -1e-14 ) q = 0;
      if (p<0 || q<0) 
      {
        return 0.0;
        std::cout << p << " p | q " << q << "\n";
        std::cout << K << std::endl; 
        std::cout << "something went wrong in Eigenvalues for beta!" << std::endl;
        assert(false);
        abort();
      }
      q = sqrt(q);
      double l_max = p + q;
      double l_min = p - q;
      
      return SQR(l_max)/l_min;
    }

    // --factorBeta
    double factorBeta(const double intelFactor, 
                      const double faceVol,
                      const FluxRangeType& enK, 
                      const FluxRangeType& nbK,
                      double& wspeedL ) const
    {
      //double minEn = 1e308;
      //double maxEn = -1e308;
      //double minNb = 1e308;
      //double maxNb = -1e308;
      double minK = 1e308;
      double maxK = -1e308;

      FluxRangeType K;

      for(int j=0; j<dimDomain; ++j) 
      {
        for (int i=0;i<dimDomain;++i) 
          K[i][j] = 0.5*(enK[i][j]+nbK[i][j]);

        minK = std::min( K[j][j] , minK); 
        maxK = std::max( K[j][j] , maxK); 

        /*
        minEn = std::min( enK[j][j] , minEn); 
        maxEn = std::max( enK[j][j] , maxEn); 
        
        minNb = std::min( nbK[j][j] , minNb); 
        maxNb = std::max( nbK[j][j] , maxNb); 
        */
      }

      // store local diffusion time step 
      wspeedL = std::max( std::abs(maxK) , std::abs(minK) );

      //double betEn = SQR(maxEn) / minEn;
      //double betNb = SQR(maxNb) / minNb;
      //double betS  = SQR(maxK) / minK;
      //double betS  = std::abs(maxK) / minK;
       
      //compBetaK(enK,betEn);
      //compBetaK(nbK,betNb);
      double betS = compBetaK(K);
      
      //double betS  = SQR(maxK) / minK;
      //double betS  = 1./ SQR(maxK);// / minK;
      //double betS  = std::abs(maxK) / minK;
       
      /*
      //std::cout << betS << " betS factor vorher \n";
      if (&enK != &nbK) 
      {
        double jumpK = tanh(std::abs(betEn-betNb));
        //std::cout << jumpK << " jump \n";
        betS = betS * jumpK + betS * (1.-jumpK);
      }
      */

      //std::cout << betS << " betS factor nachher \n";
      return (beta_ * intelFactor * betS * faceVol);
    }

    template <class QuadratureImp> 
    double applyLocalNeighbor(const IntersectionType & nit, 
                              EntityType & en, 
                              EntityType & nb,
                              VolumeQuadratureType & volQuad,
                              const QuadratureImp & faceQuadInner, 
                              const QuadratureImp & faceQuadOuter, 
                              const BaseFunctionSetType & bsetEn, 
                              LocalMatrixType & matrixEn,
                              SingleLFType& singleRhs,
                              double& wspeedS 
#ifdef DG_DOUBLE_FEATURE
                              , const bool interior 
#endif    
                              ) const
    {
      const int numDofs = bsetEn.numBaseFunctions();

      // make neighbor known to model caller 
      caller_.setNeighbor(nb);

      ////////////////////////////////////////////////////////////
      RangeType resultLeft(0.0);
      RangeType resultRight(0.0);

      RangeType phi_j;
      RangeType phiNeigh_j;
      RangeType phiEn;
      RangeType phiNeigh;

      // reuse cache mem 
      assert( psi_.size() > 0 );
      assert( coeffPsi_.size() > 0 );
      JacobianRangeType& normEn = psi_[0];
      JacobianRangeType& normNb = coeffPsi_[0];

      // create matrix handles for neighbor 
      LocalMatrixType matrixNb = matrixObj_.localMatrix( en, nb );

#ifdef DG_DOUBLE_FEATURE
      // create matrix handles for neighbor (when called with ghost do nothing)
      LocalMatrixType enMatrix = matrixObj_.localMatrix( nb, (interior) ? en : nb ); 

      // create matrix handles for neighbor 
      LocalMatrixType nbMatrix = matrixObj_.localMatrix( nb, nb ); 
      
      // set matrix to id matrix 
      if( ! interior ) 
      {
        for(int k=0; k<numDofs; ++k) 
        {
          nbMatrix.set( k, k, 1);
        }
      }
#else 
      bool useInterior = false;
#endif
      // get base function set 
      const BaseFunctionSetType bsetNeigh = spc_.baseFunctionSet(nb);

      typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
      const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
      const BaseFunctionSetType nbSet = gradientSpace_.baseFunctionSet( nb );

      typedef FieldMatrix<double, massSize , massSize > MassMatrixType; 
      typedef FieldVector<double, massSize > MassVectorType; 

      // loop over all quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // cacluate outer normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol; 

        // make sure we have the same factors 
        assert( std::abs(faceQuadInner.weight(l) - faceQuadOuter.weight(l)) < 1e-10);
        // integration element factor 
        const double intelFactor = faceQuadInner.weight(l);
        
        const double intel = faceVol * intelFactor; 
#ifdef DG_DOUBLE_FEATURE
        // use opposite signs here
        const double outerIntel = -faceVol * faceQuadOuter.weight(l); 
        // intel switching between bilinear from B_+ and B_-  
        const double outerBilinIntel = (bilinearPlus_) ? outerIntel : -outerIntel;

        // we alwas stay on the positive side 
        const RangeFieldType C_12 = 0.5;
#else 
#endif
        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        if(problem_.hasCoefficient())
        {
          // call anayltical flux of discrete model 
          caller_.evaluateCoefficientFace(nit,
              faceQuadInner,faceQuadOuter,l,coeffEn_,coeffNb_);

          for(int i=0; i<dimRange; ++i)
          {
            normEn[i] = 0.0;
            normNb[i] = 0.0;
            coeffEn_.umv(unitNormal,normEn[i]);
            coeffNb_.umv(unitNormal,normNb[i]);
          }
        }
        else 
        {
          // set to unit matrix 
          setToUnitMatrix(coeffEn_, coeffNb_);

          for(int i=0; i<dimRange; ++i)
          {
            for(int j=0; j<dimDomain; ++j)
            {
              normEn[i][j] = unitNormal[j];
              normNb[i][j] = unitNormal[j];
            }
          }
        }

        double ldt = 0.0;
        // overall beta factor 
        const double facBeta = factorBeta(intelFactor,faceVol,coeffEn_,coeffNb_, ldt);

        wspeedS += ldt * faceQuadInner.weight(l) ;

        // C_12 switch 
        const RangeFieldType C_12 = ((unitNormal * upwind_) < 0) ? -0.5 : 0.5;
        // set useInterior to save comp time 
        useInterior = ( C_12 > 0 );
        
        // cache base functions evaluations
        // leads to major speedup
        for(int k=0; k<numDofs; ++k)
        { 
          // eval base functions 
          bsetEn.evaluate(k,faceQuadInner[l], phi_[k]);
          // eval gradient for en 
          tau_[k] = bsetEn.evaluateGradientSingle(k, en, faceQuadInner[l] , normEn);  

          // neighbor stuff 
          bsetNeigh.evaluate(k,faceQuadOuter[l], phiNeigh_[k] );      
          // eval gradient for nb 
          tauNeigh_[k] = bsetNeigh.evaluateGradientSingle(k, nb, faceQuadOuter[l] , normNb);      
        }
               
        // this terms dissapear if Babuska-Zlamal is used 
        {
          for(int k=0; k<numDofs; ++k)
          {
            for(int j=0; j<numDofs; ++j)
            {
              // view from inner entity en 
              // v^+ * (grad w^+  + grad w^-)
              {
                numericalFlux2(phi_[k] , tau_[j] , tauNeigh_[j] , resultLeft, resultRight);

                RangeFieldType valLeft = resultLeft[0];
                valLeft *= -intel;

                matrixEn.add( k , j , valLeft );

                RangeFieldType valRight = resultRight[0];
                valRight *= -intel;

                matrixNb.add( k , j , valRight );
              }

              // view from inner entity en 
              // grad v^+ * ( w^+  - w^-)
              {
                numericalFlux(tau_[k] , phi_[j] , phiNeigh_[j] , resultLeft, resultRight);

                RangeFieldType valLeft = resultLeft[0];
                valLeft *= bilinIntel;

                matrixEn.add( k , j , valLeft );

                RangeFieldType valRight = resultRight;
                valRight *= bilinIntel;

                matrixNb.add( k , j , valRight );
              }
#ifdef DG_DOUBLE_FEATURE 
              // this part should only be calculated if neighboring
              // entity has partition type interior 
              if( interior ) 
              {
                // view from outer entity nb 
                // v^+ * (grad w^+  + grad w^-)
                {
                  numericalFlux2(phiNeigh_[k] , tauNeigh_[j] , tau_[j] , resultLeft, resultRight);

                  RangeFieldType valLeft = resultLeft[0];
                  valLeft *= -outerIntel;

                  nbMatrix.add( k , j , valLeft );

                  RangeFieldType valRight = resultRight[0];
                  valRight *= -outerIntel;

                  enMatrix.add( k , j , valRight );
                }

                // view from outer entity nb 
                // v^+ * (grad w^+  + grad w^-)
                {
                  numericalFlux(tauNeigh_[k] , phiNeigh_[j] , phi_[j] , resultLeft, resultRight);

                  RangeFieldType valLeft = resultLeft;
                  valLeft *= outerBilinIntel;

                  nbMatrix.add( k , j , valLeft );

                  RangeFieldType valRight = resultRight;
                  valRight *= outerBilinIntel;

                  enMatrix.add( k , j , valRight );
                }
              }
#endif
            } // end for 
          } // end for 
        } // end }  

        if( betaNotZero_ )
        {
          for(int k=0; k<numDofs; ++k)
          {
            for(int j=0; j<numDofs; ++j)
            {
              // phi_j * phi_k on entity 
              phi_j    = phi_[j] * phi_[k]; 
              // product with nb 
              phiNeigh = phiNeigh_[j] * phi_[k]; //bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phi_[k] );      
              
              // phi_j * phi_k on neighbour  
              phiNeigh_j = phiNeigh_[j] * phiNeigh_[k];//bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phiNeigh_[k] );      
              // product with nb 
              phiEn = phi_[j] * phiNeigh_[k]; //bsetEn.evaluateSingle(j,faceQuadInner,l, phiNeigh_[k]); 

              // view from inner entity en 
              {
                numericalFluxStab(phi_j, phiNeigh , resultLeft, resultRight);

                RangeFieldType valLeft = facBeta;
                valLeft  *= resultLeft;

                matrixEn.add( k , j , valLeft );

                RangeFieldType valRight = facBeta;
                valRight *= resultRight;

                matrixNb.add( k , j , valRight );
              }
              
#ifdef DG_DOUBLE_FEATURE
              // view from outer entity nb
              {
                numericalFluxStab(phiNeigh_j, phiEn , resultLeft, resultRight);

                RangeFieldType valLeft = facBeta;
                valLeft  *= resultLeft;

                nbMatrix.add( k , j , valLeft );

                RangeFieldType valRight = facBeta;
                valRight *= resultRight;

                enMatrix.add( k , j , valRight );
              }
#endif
            }
          }
        }
      } // end loop quadrature points 

      if( timeDependent_ && theta_ < 1.0 ) 
      {
        assert( uh_ );
        const SingleLFType nbLf = uh_->localFunction( nb );
        multLocal( matrixNb, singleRhs, nbLf );
      }

      return nb.geometry().volume();
    } // end applyLocalNeighbor 

    template <class BaseFunctionSet, 
              class LocalStorageType, 
              class MassMatrixType>
    void getMassMatrix(const GeometryType& geo,
                       VolumeQuadratureType& volQuad,
                       BaseFunctionSet& set,
                       const int numBase,
                       LocalStorageType& tmp,
                       MassMatrixType& massMatrix) const
    {
      const int volNop = volQuad.nop();
      for(int qp=0; qp<volNop; ++qp) 
      {
        // calculate integration weight 
        const double intel = volQuad.weight(qp)
           * geo.integrationElement(volQuad.point(qp));

        for(int m=0; m<numBase; ++m)
        {  
          // eval base functions 
          set.evaluate(m, volQuad[qp], tmp[m] );
        }

        for(int m=0; m<numBase; ++m)
        {
          {
            double val = intel * (tmp[m] * tmp[m]);
            massMatrix[m][m] += val;
          }

          
          for(int k=m+1; k<numBase; ++k) 
          {
            double val = intel * (tmp[m] * tmp[k]);
            massMatrix[k][m] += val;
            massMatrix[m][k] += val;
          }
        }
      }
    }

  private:  
    void numericalFlux(const RangeFieldType & grad, 
                       const RangeType & phiLeft,
                       const RangeType & phiRight, 
                       RangeType & resultLeft,
                       RangeType & resultRight) const
    {
      resultLeft  = grad;
      resultLeft *= 0.5;
      
      resultRight  = grad;
      // need negative value of phiRight 
      resultRight *= -0.5; 

      resultLeft  *= phiLeft;
      resultRight *= phiRight;
    }
                       
    void numericalFlux2(const RangeType & phi,
                        const RangeFieldType & gradLeft, 
                        const RangeFieldType & gradRight,
                        RangeType & resultLeft,
                        RangeType & resultRight) const
    {
      resultLeft  = gradLeft;
      resultLeft *= 0.5 * phi[0];
      
      resultRight  = gradRight;
      resultRight *= 0.5 * phi[0]; 
    }
                       
    void numericalFluxStab(const RangeType & phiLeft,
                           const RangeType & phiRight,
                           RangeType & resultLeft,
                           RangeType & resultRight) const
    {
      resultLeft  = phiLeft; 
      resultRight = phiRight; 
      resultRight *= -1.0;
    }

    void numericalFlux_C12(const RangeType & phi,
                        const RangeFieldType & gradLeft,
                        const RangeFieldType & gradRight,
                        RangeType & resultLeft,
                        RangeType & resultRight) const
    {
      resultLeft  =  gradLeft;
      resultRight = -gradRight;

      resultLeft  *= phi;
      resultRight *= phi;
    }

    void numericalFlux2_C12(const RangeFieldType & grad,
                            const RangeType & phiLeft,
                            const RangeType & phiRight,
                            RangeType & resultLeft,
                            RangeType & resultRight) const
    {
      resultLeft  =  phiLeft;
      resultRight = -phiRight;

      resultLeft  *= grad;
      resultRight *= grad;
    }
    
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    //! empty constructor not defined 
    DGPrimalOperatorImpl();
    //! copy constructor not defined 
    DGPrimalOperatorImpl(const DGPrimalOperatorImpl&);

  protected:
    mutable DiscreteModelCallerType caller_;

    DiscreteModelType& problem_; 
         
    
    mutable ArgumentType* arg_;
    mutable DestinationType* rhs_;
    mutable const DestinationType* uh_;

    const DiscreteFunctionSpaceType& spc_;
    GridPartType & gridPart_;
    const DiscreteGradientSpaceType gradientSpace_;

    const LocalIdSetType & localIdSet_;
    const GridWidthType& gridWidth_;
    
    const int volumeQuadOrd_;
    const int faceQuadOrd_;

    mutable MatrixObjectType matrixObj_;

    // return type of analyticalFlux 
    mutable FluxRangeType coeffEn_;
    mutable FluxRangeType coeffNb_;
    // caches for base function evaluation 
    mutable MutableArray<RangeFieldType> tau_;
    mutable MutableArray<RangeFieldType> tauNeigh_;
    mutable MutableArray<RangeType> phi_;
    mutable MutableArray<RangeType> phiNeigh_;
    mutable MutableArray<JacobianRangeType> psi_;
    mutable MutableArray<JacobianRangeType> coeffPsi_;

    mutable MutableArray<RangeType> eta_;
    mutable MutableArray<RangeType> etaNeigh_;

    mutable MutableArray<GradRangeType> rRets_;
    mutable MutableArray<GradRangeType> rRetsCoeff_;

    DomainType upwind_;

    mutable JacobianRangeType psitmp_;

    mutable bool matrixAssembled_;
    double betaFactor_;
    mutable double globalBeta_;
    mutable double beta_;

    // if true B_+ is used otherwise B_-
    bool bilinearPlus_;
    double power_;
    bool betaNotZero_;
    mutable double factorFaces_;
    mutable double tau_1_;
    mutable double theta_;
    mutable double dtMin_;
    const double minLimit_;
    bool timeDependent_;
    int sequence_ ;
  };

  template <class DiscreteModelImp, 
            class PreviousPassImp, 
            class MatrixObjectTraits, 
            int passId = -1 >
  class DGPrimalOperator
    : public DGPrimalOperatorImpl< DiscreteModelImp, 
                                    PreviousPassImp, 
                                    MatrixObjectTraits, 
                                    passId > ,
     public OEMSolver :: PreconditionInterface                               
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef DGPrimalOperatorImpl<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > BaseType;

    typedef DGPrimalOperator<DiscreteModelImp,
            PreviousPassImp,MatrixObjectTraits, passId > ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;

    typedef PreviousPassImp PreviousPassType;

    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    typedef typename BaseType :: MatrixObjectType MatrixObjectType;
    typedef typename BaseType :: MatrixObjectType::MatrixType MatrixType;
    typedef typename BaseType :: MatrixObjectType::PreconditionMatrixType PreconditionMatrixType;
    
  protected:  
    using BaseType :: matrixObj_ ;
  public:
    //- Public methods
    /**  \brief Constructor
     \param problem Actual problem definition (see problem.hh)
     \param pass Previous pass
     \param spc Space belonging to the discrete function local to this pass
     \param paramFile parameter file to read necessary parameters, if empty 
             default parameters will be applied 
    
     \note Available methods are (chosen by parameters B_{+,-}, beta)
          - Interior Penalty : B_{+,-}: 0 , beta: > 0 (big)
          - Baumann-Oden     : B_{+,-}: 1 , beta: = 0       
          - NIPG             : B_{+,-}: 1 , beta: > 0       
     */         
    DGPrimalOperator(DiscreteModelType& problem, 
                      PreviousPassType& pass, 
                      const DiscreteFunctionSpaceType& spc,
                      const std::string paramFile = "")
      : BaseType(problem, pass, spc, paramFile)
    {
    }
  public:
    virtual ~DGPrimalOperator() {}

    //! return refernence to system matrix, used by Solvers
    const MatrixObjectType & systemMatrix () const { return matrixObj_; }
    
    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType & preconditionMatrix () const { 
      return matrixObj_.preconditionMatrix(); 
    }

    //! returns true if preconditioning matrix has been build 
    bool hasPreconditionMatrix() const  { 
      return matrixObj_.hasPreconditionMatrix(); 
    }
  };
#undef DG_DOUBLE_FEATURE  
} // end namespace Dune
#endif
