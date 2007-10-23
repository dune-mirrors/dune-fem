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

#include <dune/fem/misc/timeprovider.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/function/common/dfcommunication.hh>
#include <dune/fem/space/common/communicationmanager.hh>
#include <dune/fem/space/common/arrays.hh>

#include <dune/fem/misc/gridwidth.hh>

#include "dgmatrixsetup.hh"

// double feature only works in serial runs 
//#if HAVE_MPI == 0
#define DG_DOUBLE_FEATURE 
//#endif

namespace Dune {
/*! @ingroup EllipticOperator
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) + &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/
  ////////////////////////////////////////////////////////////
  //
  //  --DGPrimalOperator 
  //
  ////////////////////////////////////////////////////////////
  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class GradientPassImp, 
            class PreviousPassImp, class MatrixObjectImp>
  class DGPrimalOperator 
    : public LocalPass<DiscreteModelImp, PreviousPassImp> 
    , public OEMSolver::PreconditionInterface
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    typedef DGPrimalOperator<DiscreteModelImp,GradientPassImp,
            PreviousPassImp,MatrixObjectImp> ThisType;

    typedef GradientPassImp GradientPassType; 
    typedef typename GradientPassType :: DiscreteModelType
      GradientDiscreteModelType;

    //! space of gradients of function 
    typedef typename GradientPassType :: 
      DiscreteFunctionSpaceType DiscreteGradientSpaceType;

    typedef typename GradientPassType :: DestinationType GradDestinationType;
    typedef typename GradientPassType :: DiscreteModelCallerType GradientModelCallerType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
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
        
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType :: GridType :: Traits :: LocalIdSet LocalIdSetType; 
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef FieldVector<double,DomainType::dimension-1> FaceDomainType;

    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef EllipticDiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    enum { dim = GridType :: dimension };
    
    typedef typename GridType :: ctype ctype;
    typedef FieldMatrix<ctype,dim,dim> JacobianInverseType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType::polynomialOrder };

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
    
    //! type of underlying matrix implementation 
    typedef MatrixObjectImp MatrixObjectType; 

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

        for (int j = 0; j < numDofs_; ++j) 
        {
          singleRhs_[j] += bsetEn_.evaluateSingle(j, volQuad, l, rhsval_ );
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
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param paramFile parameter file to read necessary parameters, if empty 
    //!         default parameters will be applied 
    //!
    //!  NOTE: possible parameters are
    //!     - beta if beta = 0 then Baumann-Oden is chosen, otherwise beta > 0, default is 0
    //!     - B_{+,-} choose between B_+ and B_-, 1 == B_+ | 0 == B_- ,default is B_+
    //!     - Babuska-Zlamal 1 means we take Babuska-Zlamal method, 0 not, default is 0 
    //!       if Babuska-Zlamal is chosen, beta > 0 is needed
    //!         
    DGPrimalOperator(DiscreteModelType& problem, 
                GradientPassType & gradPass,
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                const std::string paramFile = "")
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      localIdSet_(gridPart_.grid().localIdSet()),
      gridWidth_ ( GridWidthProviderType :: getObject( &spc_.grid())),
      time_(0),
      volumeQuadOrd_(2* spc_.order() ),
      faceQuadOrd_(2*spc_.order() + 1),
      matrixObj_(spc_,spc_, paramFile ),
      upwind_(M_PI),
      coeffEn_(1.0),
      coeffNb_(1.0),
      matrixAssembled_(false),
      betaFactor_(0.0),
      globalBeta_(1.0),
      beta_(0.0),
      bilinearPlus_(true),
      power_( 2*spc_.order() ),
      notBabuskaZlamal_(true),
      compactLDG_(false),
      betaNotZero_(false)
    {
      if( dim > 1 ) upwind_[1] = M_LN2;
      if( dim > 2 ) upwind_[2] = M_E;

      if( ! (spc_.order() > 0))
      {
        std::cerr << "ERROR: DG Primal operator only working for spaces with polynomial order > 0! \n";
        assert(false);
        abort();
      }
      
      // if parameter file is not empty read parameter 
      if(paramFile != "")
      {
        const bool output = (gridPart_.grid().comm().rank() == 0);
        readParameter(paramFile,"beta",betaFactor_, output);
        int bplus = 1;
        readParameter(paramFile,"B_{+,-}",bplus, output);
        assert( (bplus == 0) || (bplus == 1) ); 
        bilinearPlus_ = (bplus == 0) ? false : true; 

        int zlamal = 0;
        readParameter(paramFile,"Babuska-Zlamal",zlamal, output);
        notBabuskaZlamal_ = (zlamal == 1) ? false : true;
      }

      betaNotZero_ = (std::abs(betaFactor_) > 0.0);

      /*
      if( !betaNotZero_ && !notBabuskaZlamal_)
      {
        std::cerr << "ERROR: beta == 0.0 and Babuska-Zlamal == 1 !";
        std::cerr << " Choose either beta > 0 or Babuska-Zlamal = 1" << std::endl;
        assert(false);
        exit(1);
      }
      */

      if( ! betaNotZero_ && notBabuskaZlamal_ )
      {
        std::cout << "DGPrimalOperator: using Baumann-Oden method!\n"; 
        if(spc_.order() < 2)
        {
          std::cerr << "WARNING: Baumann-Oden method only working for polynomial degree >= 2! \n";
          assert(false);
          exit(1);
        }
      }
      else 
      {
        if( notBabuskaZlamal_ )
        {
          if(bilinearPlus_ )
          {
            std::cout << "DGPrimalOperator: using NIPG method, beta = " << betaFactor_ << " !\n"; 
          }
          else 
          {
            std::cout << "DGPrimalOperator: using Interior Penalty method, beta = " << betaFactor_ << " !\n"; 
          }
        }
        else 
        {
          if(bilinearPlus_ )
          {
            std::cout << "DGPrimalOperator: using Babuska-Zlamal method, beta = " << betaFactor_ << " !\n"; 
          }
          else 
          {
            notBabuskaZlamal_ = true;
            compactLDG_ = true; 
            if( ! betaNotZero_ )
            {
              DUNE_THROW(InvalidStateException,"Beta > 0 for Compact LDG not provided!");
            }
            std::cout << "DGPrimalOperator: using Compact LDG method, beta = " << betaFactor_ << " !\n"; 
          }
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

    //! Destructor
    virtual ~DGPrimalOperator() 
    {
    }

    //! Stores the time provider passed by the base class in order to have
    //! access to the global time
    virtual void processTimeProvider(TimeProvider* time) {
      time_ = time;
    }

    //! Estimate for the timestep size, return value is 0 
    double timeStepEstimate() const {
      return 0.0;
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, DestinationType & rhs)
    {
      // prepare 
      prepare( arg, rhs );

      // if grid has changed, then matrix structure is re-build
      {
        // object for creating the stencil 
        ElementAndNeighbors stencil;
        matrixObj_.reserve(stencil);
      }

      // clear matrix 
      matrixObj_.clear();
      
      // clear right hand side 
      rhs.clear();

      // build matrix and rhs 
      //Timer timer; 
      this->compute( arg, rhs );
      matrixAssembled_ = true;
      //std::cout << "Setup of Matrix took " << timer.elapsed() << " sec.\n";

      // create pre-condition matrix if activated 
      matrixObj_.createPreconditionMatrix();
      
      // finalize 
      finalize( arg, rhs );
    }

  public:  
    //! do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer());
    }
    
    //! do matrix vector multiplication, used by OEM-Solver and DuneODE Solvers  
    // --multOEM
    void multOEM(const double * arg, double * dest) const
    {
      matrixObj_.multOEM( arg, dest );
    }

  public:
    //! return refernence to system matrix, used by Solvers
    const MatrixObjectType & systemMatrix () const { return matrixObj_; }
    
    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType & preconditionMatrix () const { 
      return matrixObj_.pcMatrix(); 
    }

    //! returns true if preconditioning matrix has been build 
    bool hasPreconditionMatrix() const  { 
      return matrixObj_.hasPcMatrix(); 
    }
           
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      caller_.setArgument(*arg_);

      // calculate beta = O(1/h)
      globalBeta_ = betaFactor_/gridWidth_.gridWidth();

      if (time_) {
        caller_.setTime(time_->time());
      }
      else {
        caller_.setTime(0.0);
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      caller_.finalize();
    }

  public:
    //! return reference to function space 
    const DiscreteFunctionSpaceType & space() const 
    {
      return spc_;
    }

  public:
    //! set all data belonging to this entity to zero 
    void restrictLocal(const EntityType& father, const EntityType & son, bool firstCall ) const
    {
    }

    //! set all data belonging to this entity to zero 
    void prolongLocal(const EntityType& father, const EntityType & son, bool ) const
    {
    }

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

      // set default value to fMat 
      coeffEn_ = 1.0;
      coeffNb_ = 1.0;
      
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
          bsetEn.jacobian( k, volQuad, l, psitmp_ );
  
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
      
    //! apply operator on entity 
    void applyLocal(EntityType& en) const
    {
      // only build Matrix in interior 
      assert( en.partitionType() == InteriorEntity );
      
      // get local element matrix 
      LocalMatrixType matrixEn = matrixObj_.localMatrix(en,en); 

      // make entities known in callers
      caller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      const GeometryType & geo = en.geometry();

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );

      // local function for right hand side 
      SingleLFType singleRhs = dest_->localFunction(en); //rhs
      
      // resize caches 
      resizeCaches(numDofs);

      double betaEst = 0.0;
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////
      if(problem_.hasCoefficient() && problem_.hasRHS() )
      {
        CoefficientCaller<DiscreteModelCallerType,true,true> coeffCaller( singleRhs ); 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else if( problem_.hasCoefficient() )
      {
        CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else if ( problem_.hasRHS() )
      {
        CoefficientCaller<DiscreteModelCallerType,false,true> coeffCaller( singleRhs ); 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }
      else 
      {
        CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
        betaEst = volumetricPart(en,geo,volQuad,coeffCaller,bsetEn,numDofs,matrixEn);
      }

      beta_ = betaEst * globalBeta_;

      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // if neighbor exists 
        if (nit.neighbor()) 
        {
          EntityPointerType neighEp = nit.outside();
          EntityType& nb = *neighEp;

#ifdef DG_DOUBLE_FEATURE
          // get partition type 
          const bool nonInterior = 
            ( nb.partitionType() != InteriorEntity );
          // only once per intersection or when outside is not interior 
          if( (localIdSet_.id(en) < localIdSet_.id(nb)) 
              || nonInterior
            )
#endif
          {
            // type of TwistUtility 
            typedef TwistUtility<GridType> TwistUtilityType;
            // check conformity 
            if( TwistUtilityType::conforming(gridPart_.grid(),nit) ) 
            {
              FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                               FaceQuadratureType::INSIDE);

        
              FaceQuadratureType faceQuadOuter(gridPart_, nit, faceQuadOrd_,
                                               FaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              applyLocalNeighbor(nit,en,nb,volQuad,
                    faceQuadInner,faceQuadOuter, 
                    bsetEn,matrixEn
#ifdef DG_DOUBLE_FEATURE
                    , ! nonInterior 
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
                nonConformingFaceQuadInner(gridPart_, nit, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::INSIDE);
          
              NonConformingFaceQuadratureType 
                nonConformingFaceQuadOuter(gridPart_,nit, faceQuadOrd_,
                                           NonConformingFaceQuadratureType::OUTSIDE);

              // apply neighbor part 
              applyLocalNeighbor(nit,en,nb,volQuad,
                    nonConformingFaceQuadInner,
                    nonConformingFaceQuadOuter, 
                    bsetEn,matrixEn
#ifdef DG_DOUBLE_FEATURE
                    , ! nonInterior 
#endif
                    );
            }
          }
        } // end if neighbor 

        // if intersection with boundary 
        if (nit.boundary()) 
        { 
          // create quadrature 
          FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                           FaceQuadratureType::INSIDE);

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

            // overall beta factor 
            const double facBeta = factorBeta(intelFactor,faceVol);

            // get boundary value 
            RangeType boundaryValue(0.0);

            // call boundary value function 
            BoundaryIdentifierType bndType = 
              caller_.boundaryValue(nit,faceQuadInner,l,boundaryValue);

            // only Dirichlet and Neumann Boundary supported right now 
            assert( bndType.isDirichletType() || bndType.isNeumannType() );

            ///////////////////////////////
            //  evaluate coefficients 
            ///////////////////////////////
            assert( psi_.size() > 0 );
            JacobianRangeType& norm = psi_[0];
            if(problem_.hasCoefficient())
            {
              // evaluate coefficient on boundary
              caller_.evaluateCoefficientBoundary(nit,faceQuadInner,l,coeffEn_);
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

            // cache base functions evaluations
            for(int k=0; k<numDofs; ++k)
            { 
              // evaluate normal * grad phi 
              tau_[k] = bsetEn.evaluateGradientSingle(k,en,faceQuadInner,l, norm);  
              // evaluate phi 
              bsetEn.evaluate(k,faceQuadInner,l, phi_[k]);
            }
               
            // if not Babuska-Zlamal method, add boundary terms 
            if( notBabuskaZlamal_ )
            {
              if( bndType.isDirichletNonZero())
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
              if( bndType.isDirichletType() )
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
            if(bndType.isNeumannNonZero())
            {
              // fill matrix entries 
              for(int k=0; k<numDofs; ++k)
              {  
                // only valid for dim range = 1
                RangeFieldType rhsVal = boundaryValue * phi_[k];

                rhsVal *= intelFactor;
                //rhsVal *= intel;
                singleRhs[k] += rhsVal;
              }
            }

            if( betaNotZero_ )
            {
              // stabilization 
              if( bndType.isDirichletType())
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
              }
                
              // dirichlet boundary values for u 
              if(bndType.isDirichletNonZero())
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
        } // end if boundary

      } // end intersection iterator 

      // resort corresponding matrix rows for ascending numbering 
      matrixEn.resort(); 
    } // end apply local 

    double factorBeta(const double intelFactor, const double faceVol) const 
    {
      if(notBabuskaZlamal_)
      {
        return (beta_ * intelFactor);
      }
      else 
      {
        return (beta_ * intelFactor * faceVol * pow(faceVol , -power_ )); 
      }
    }

    template <class QuadratureImp> 
    void applyLocalNeighbor(IntersectionIteratorType & nit, 
                            EntityType & en, 
                            EntityType & nb,
                            VolumeQuadratureType & volQuad,
                            const QuadratureImp & faceQuadInner, 
                            const QuadratureImp & faceQuadOuter, 
                            const BaseFunctionSetType & bsetEn, 
                            LocalMatrixType & matrixEn
#ifdef DG_DOUBLE_FEATURE
                            , const bool interior 
#endif
                            ) const
    {
      const int numDofs = bsetEn.numBaseFunctions();

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
      // create matrix handles for neighbor 
      LocalMatrixType enMatrix = matrixObj_.localMatrix( nb, en ); 
      
      // create matrix handles for neighbor 
      LocalMatrixType nbMatrix = matrixObj_.localMatrix( nb, nb ); 
#endif
      // get base function set 
      const BaseFunctionSetType bsetNeigh = spc_.baseFunctionSet(nb);
     
      // loop over all quadrature points 
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // cacluate outer normal 
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol; 

        // C_12 stabilization factor 
        const RangeFieldType C_12 = (unitNormal * upwind_ < 0) ? -0.5 : 0.5;

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
#endif
        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        // overall beta factor 
        const double facBeta = factorBeta(intelFactor,faceVol);

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
          for(int i=0; i<dimRange; ++i)
          {
            for(int j=0; j<dimDomain; ++j)
            {
              normEn[i][j] = unitNormal[j];
              normNb[i][j] = unitNormal[j];
            }
          }
        }
               
        // cache base functions evaluations
        // leads to major speedup
        for(int k=0; k<numDofs; ++k)
        { 
          // eval base functions 
          bsetEn.evaluate(k,faceQuadInner,l, phi_[k]);
          // eval gradient for en 
          tau_[k] = bsetEn.evaluateGradientSingle(k, en, faceQuadInner, l, normEn);  

          // neighbor stuff 
          bsetNeigh.evaluate(k,faceQuadOuter,l, phiNeigh_[k] );      
          // eval gradient for nb 
          tauNeigh_[k] = bsetNeigh.evaluateGradientSingle(k, nb, faceQuadOuter, l, normNb);      
        }
               
        // this terms dissapear if Babuska-Zlamal is used 
        if(notBabuskaZlamal_)
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
              // v^+ * (grad w^+  + grad w^-)
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

              ////////////////////////////////////
              //  C12 stabilization 
              ///////////////////////////////////

              if(compactLDG_)
              {
                // view from inner entity en 
                {
                  numericalFlux_C12(phi_[k], tau_[j] , tauNeigh_[j] , resultLeft, resultRight);

                  RangeFieldType valLeft = (C_12 * resultLeft[0]);
                  valLeft *= bilinIntel;

                  matrixEn.add( k , j , valLeft );

                  RangeFieldType valRight = (C_12 * resultRight[0]);
                  valRight *= bilinIntel;

                  matrixNb.add( k , j , valRight );
                }

                // view from inner entity en 
                {
                  numericalFlux2_C12(tau_[k] , phi_[j] , phiNeigh_[j] , resultLeft, resultRight);

                  RangeFieldType valLeft = (resultLeft[0] * C_12);
                  valLeft *= bilinIntel;

                  matrixEn.add( k , j , valLeft );

                  RangeFieldType valRight = (resultRight[0] * C_12);
                  valRight *= bilinIntel;

                  matrixNb.add( k , j , valRight );
                }

#ifdef DG_DOUBLE_FEATURE 
                // this part should only be calculated if neighboring
                // entity has partition type interior 
                if( interior ) 
                {
                  // view from inner entity en 
                  {
                    numericalFlux_C12(phiNeigh_[k], tauNeigh_[j] , tau_[j] , resultLeft, resultRight);

                    RangeFieldType valLeft = (C_12 * resultLeft[0]);
                    valLeft *= bilinIntel;

                    nbMatrix.add( k , j , valLeft );

                    RangeFieldType valRight = (C_12 * resultRight[0]);
                    valRight *= bilinIntel;

                    enMatrix.add( k , j , valRight );
                  }

                  // view from inner entity en 
                  {
                    numericalFlux2_C12(tauNeigh_[k] , phiNeigh_[j] , phi_[j] , resultLeft, resultRight);

                    RangeFieldType valLeft = (resultLeft[0] * C_12);
                    valLeft *= bilinIntel;

                    nbMatrix.add( k , j , valLeft );

                    RangeFieldType valRight = (resultRight[0] * C_12);
                    valRight *= bilinIntel;

                    enMatrix.add( k , j , valRight );
                  }
                }
#endif
              } // compact LDG 
            } // end for 
          } // end for 
        } // end notBabuskaZlamal

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
    DGPrimalOperator();
    //! copy constructor not defined 
    DGPrimalOperator(const DGPrimalOperator&);

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const LocalIdSetType & localIdSet_;
    const GridWidthType& gridWidth_;
    
    // time provider 
    TimeProvider* time_;

    const int volumeQuadOrd_;
    const int faceQuadOrd_;

    mutable MatrixObjectType matrixObj_;
    DomainType upwind_;

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
    mutable JacobianRangeType psitmp_;

    mutable bool matrixAssembled_;
    double betaFactor_;
    mutable double globalBeta_;
    mutable double beta_;

    // if true B_+ is used otherwise B_-
    bool bilinearPlus_;
    double power_;
    bool notBabuskaZlamal_;
    bool compactLDG_;
    bool betaNotZero_;
  };
#undef DG_DOUBLE_FEATURE  
} // end namespace Dune
#endif
