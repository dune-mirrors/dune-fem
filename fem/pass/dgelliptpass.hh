#ifndef DUNE_ELLIPTPASS_HH
#define DUNE_ELLIPTPASS_HH

#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/modelcaller.hh>

// * needs to move
#include <dune/fem/misc/timeutility.hh>

#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/utility/twistutility.hh>
#include <dune/fem/misc/boundaryidentifier.hh>

#include <dune/common/typetraits.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/discretefunction/common/dfcommunication.hh>

#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune {

  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGElliptGradientPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    typedef LocalDGElliptGradientPass<DiscreteModelImp,PreviousPassImp> ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd =DiscreteFunctionSpaceType::polOrd};

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param quadOrd0 defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param quadOrd1 defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGElliptGradientPass(DiscreteModelType& problem, 
                    PreviousPassType& pass, 
                    DiscreteFunctionSpaceType& spc) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      spc_(spc)
    {
    }

    //! don't allocate memory here 
    virtual void allocateLocalMemory() {}
   
    //! return reference to caller 
    DiscreteModelCallerType & caller () { return caller_; }

    //! return previous pass of this pass 
    PreviousPassType & previousPass() { return this->previousPass_; }

    //! return problem for real fe pass 
    DiscreteModelType & problem () { return problem_; }

    const DiscreteFunctionSpaceType & space () const { return spc_; }
   
    //! Destructor
    virtual ~LocalDGElliptGradientPass() {
    }

    void applyLocal(EntityType& en) const
    {
    }

    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      abort();
    }

    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! Some timestep size management.
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    DiscreteFunctionSpaceType& spc_;
  };

  ////////////////////////////////////////////////////////////
  //
  //  --LocalDGElliptOperator 
  //
  ////////////////////////////////////////////////////////////
  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class GradientPassImp, 
            class PreviousPassImp>
  class LocalDGElliptOperator 
    : public LocalPass<DiscreteModelImp, PreviousPassImp> 
    , public OEMSolver::PreconditionInterface
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    typedef LocalDGElliptOperator<DiscreteModelImp,GradientPassImp,PreviousPassImp> ThisType;

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
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd =DiscreteFunctionSpaceType::polOrd};

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;
    
    typedef DiscreteFunctionSpaceType SingleDiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType:: IteratorType IteratorType ;

    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradJacobianRangeType;
    typedef GradJacobianRangeType GradientJacobianRangeType;
    enum { GradDimRange = GradientRangeType :: dimension };
    
    typedef typename DiscreteModelType :: Traits :: Traits ::template
      MatrixHandler<DiscreteFunctionSpaceType,DiscreteGradientSpaceType> ::
      MatrixHandlerType MatrixHandlerType; 

    typedef typename MatrixHandlerType::MatrixAddHandleType MatrixAddHandleType;
    typedef typename MatrixHandlerType::MatrixType MatrixType;
    typedef typename MatrixHandlerType::PreconditionMatrixType PreconditionMatrixType;
    
    typedef typename DiscreteModelType :: BoundaryIdentifierType BoundaryIdentifierType;    

    typedef CommunicationManager<DiscreteFunctionSpaceType> SingleCommunicationManagerType; 
    typedef CommunicationManager<DiscreteGradientSpaceType> GradCommunicationManagerType; 
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param quadOrd0 defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param quadOrd1 defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGElliptOperator(DiscreteModelType& problem, 
                GradientPassType & gradPass,
                PreviousPassType& pass, 
                DiscreteFunctionSpaceType& spc,
                int volumeQuadOrd =-1,int faceQuadOrd=-1) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      gradPass_(gradPass),
      gradCaller_(gradPass_.caller()),
      gradProblem_(gradPass_.problem()),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      gradientSpace_(gradPass_.space()),
      singleCommunicate_(spc_),
      gradCommunicate_(gradientSpace_),
      gradTmp_("FEPass::gradTmp",gradientSpace_),
      gradRhs_("FEPass::gradRhs",gradientSpace_),
      massTmp_(0),
      multTmp_("FEPass::multTmp",spc_),
      diag_(0),
      dtMin_(std::numeric_limits<double>::max()),
      fMat_(0.0),
      fMatNb_(0.0),
      gradSource_(0.0),
      gradSourceNb_(0.0),
      one_(0.0),
      gradFMat_(0.0),
      gradFMatNb_(0.0),
      valEn_(0.0),
      valNeigh_(0.0),
      baseEn_(0.0),
      baseNeigh_(0.0), 
      phi_(0.0),
      source_(0.0),
      sourceNb_(0.0),
      rhs_(0.0),
      tau_(0.0),
      tauneigh_(0.0),
      phiNeigh_(0.0),
      grads_(0.0),
      time_(0),
      twistUtil_(spc.grid()),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (2*spc_.order()) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (2*spc_.order()+1) : faceQuadOrd ),
      matrixHandler_(spc_,gradientSpace_,gradProblem_.hasSource(),problem_.preconditioning()),
      matrixAssembled_(false)                                                              
    {
      assert( matrixHandler_.hasMassMatrix() == gradProblem_.hasSource() );
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );

      if(gradProblem_.hasSource())
      {
        massTmp_ = new GradDestinationType ("FEPass::massTmp",gradientSpace_);
      }

      if(matrixHandler_.hasPcMatrix())
      {
        diag_ = new DestinationType("FEPass::diag",spc_);
      }

      for(int i=0; i<GradDimRange; ++i) one_[i][i] = 1.0;

      if(problem_.hasSource())
      {
        std::cerr << "Source for FEPass not supported yet. \n";
        abort();
      }
    }
   
    //! Destructor
    virtual ~LocalDGElliptOperator() 
    {
      delete diag_;
      delete massTmp_;
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

    void buildMatrix( const ArgumentType & arg, DestinationType & rhs, bool
        preCond = false )
    {
      preCond_ = preCond;
      // resize matrices 
      matrixHandler_.resizeAndClear();
     
      gradRhs_.clear();
      rhs.clear();
      
      dest_ = &rhs;
      // build matrix and rhs 
      this->compute(arg,rhs);
      matrixAssembled_ = true;

      // create pre-condition matrix if activated 
      createPreconditionMatrix();

      double * rhsPtr = gradRhs_.leakPointer();
      if(gradProblem_.hasSource())
      {
        double * massTmpPointer = massTmp_->leakPointer();
        matrixHandler_.massMatrix().multOEM(rhsPtr,massTmpPointer);
        rhsPtr = massTmpPointer;
      } 

      double * singleRhsPtr = rhs.leakPointer();
      // adjust rhs 
      double * multTmpPointer  = multTmp_.leakPointer();
      matrixHandler_.divMatrix().multOEM(rhsPtr,multTmpPointer );

      const int singleSize = spc_.size();
      for(register int i=0; i<singleSize; ++i) 
      {
        singleRhsPtr[i] -= multTmpPointer[i];      
      }
    }

    // --gradient
    template <class FuncType, class GradType> 
    void evalGradient(const FuncType & u, GradType & grad) const
    {
      grad.clear();
      if(gradProblem_.hasSource())
      {
        assert( massTmp_ );
        double * massTmpPointer = massTmp_->leakPointer();
        matrixHandler_.gradMatrix().multOEM(u.leakPointer(),massTmpPointer);
        (*massTmp_) += gradRhs_;
        matrixHandler_.massMatrix().multOEM(massTmpPointer,grad.leakPointer());
      } 
      else 
      {
        matrixHandler_.gradMatrix().multOEM(u.leakPointer(),grad.leakPointer());
        grad += gradRhs_;
      }
    }

    // do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer());
    }
    
    // do matrix vector multiplication, used by OEM-Solver  
    void multOEM(const double * arg, double * dest) const
    {
      double * multTmpPointer = multTmp_.leakPointer();
      double * gradTmpPointer = gradTmp_.leakPointer();
      
      communicate( arg );
      matrixHandler_.gradMatrix().multOEM(arg, gradTmpPointer );

      if(gradProblem_.hasSource())
      {
        assert( massTmp_ );
        double * massTmpPointer = massTmp_->leakPointer();
        matrixHandler_.massMatrix().multOEM(gradTmpPointer,massTmpPointer);
        // use mass tmp now 
        gradTmpPointer = massTmpPointer;
        // if we have mass, communicate mass 
        //doCommunicate( *massTmp_ );
        gradCommunicate_.exchange( *massTmp_ );
      }
      else 
      {
        gradCommunicate_.exchange( gradTmp_ );
        // otherwise communicate grad 
        //doCommunicate( gradTmp_ );
      }

      matrixHandler_.divMatrix().multOEM(gradTmpPointer, multTmpPointer );
      matrixHandler_.stabMatrix().multOEM(arg,dest);
      
      const int size = spc_.size();
      for(register int i=0; i<size; ++i) 
      {
        dest[i] += multTmpPointer[i];
      }
    }

  private:   
    void communicate(const double * x) const
    {
      DestinationType dest("DGEllipt::communicate_tmp",spc_,x);
      singleCommunicate_.exchange( dest );
      //doCommunicate(dest);
    }
      
    template <class DiscreteFuncType>
    void doCommunicate(DiscreteFuncType & dest) const
    {
      // if serial run, just return 
      if(gridPart_.grid().comm().size() <= 1) return;
      
      typedef DiscreteFunctionCommunicationHandler<DiscreteFuncType> DataHandleType;
      DataHandleType dataHandle(dest);
      
      gridPart_.communicate( dataHandle, InteriorBorder_All_Interface , ForwardCommunication);
    }
    
  public:
    const ThisType & systemMatrix () const { return *this; }
    const ThisType & myMatrix () const { return *this; }
    
    const PreconditionMatrixType & preconditionMatrix () const { return matrixHandler_.pcMatrix(); }
    bool hasPreconditionMatrix() const  { return matrixHandler_.hasPcMatrix(); }
           
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      caller_.setArgument(*arg_);
      gradCaller_.setArgument(*arg_);
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
      caller_.finalize();
      gradCaller_.finalize();
    }

  public:
    const DiscreteFunctionSpaceType & space() const {
      return spc_;
    }

  public:
    void prepareGlobal(const ArgumentType& arg, DestinationType& dest) const{
      prepare(arg,dest);
    }

    void finalizeGlobal() {
      finalize(*arg_,*dest_); 
      matrixAssembled_ =false;
    }

    void applyLocal(EntityType& en) const
    {
      // local function for right hand side 
      typedef typename DestinationType :: LocalFunctionType SingleLFType; 
      SingleLFType singleRhs = dest_->localFunction(en); //rhs
      
      // local function for gradient right hand side 
      typedef typename GradDestinationType :: LocalFunctionType GradLFType; 
      GradLFType gradRhs = gradRhs_.localFunction(en); //rhs
      
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

      MatrixAddHandleType stabMatrixEn(matrixHandler_.stabMatrix(),
                                     en, spc_, en, spc_ ); 
      MatrixAddHandleType gradMatrixEn(matrixHandler_.gradMatrix(),
                                     en, gradientSpace_, en, spc_ ); 
      MatrixAddHandleType divMatrixEn (matrixHandler_.divMatrix(),
                                     en, spc_, en, gradientSpace_ ); 
      
      MatrixAddHandleType massMatrixEn (matrixHandler_.massMatrix(),
                                     en, gradientSpace_, en, gradientSpace_ ); 
      
      typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
      typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;

      //- statements
      caller_.setEntity(en);
      gradCaller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      const VolumeQuadratureType massQuad(en, 0 );   
      const GeometryType & geo = en.geometry();
      const double massVolElInv = massVolumeInv(geo,massQuad);
      
      const BaseFunctionSetType& bsetEn = spc_.getBaseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      
      const GradientBaseFunctionSetType& grdbsetEn = gradientSpace_.getBaseFunctionSet(en);
      const int gradientNumDofs = grdbsetEn.numBaseFunctions();
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////

      //set default values 
      gradSource_   = 1.0;
      gradSourceNb_ = 1.0;
      // set to id matrix 
      gradFMat_   = one_;
      gradFMatNb_ = one_;

      fMat_       = 1.0;
      fMatNb_     = 1.0;
      
      const int quadNop = volQuad.nop();
      GradJacobianRangeType helpmatr(one_);

      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc factor for bas functions 
        double intel = volQuad.weight(l)*
            geo.integrationElement(volQuad.point(l))*massVolElInv;
        
        ////////////////////////////////////
        // create rightHandSide
        ////////////////////////////////////
        {
          // set default value 
          rhs_ = 0.0;

          // eval rightHandSide function 
          // if empty, rhs stays 0.0
          caller_.rightHandSide(en, volQuad, l, rhs_ );

          // scale with quadrature weight  
          rhs_ *= volQuad.weight(l);

          for (int j = 0; j < numDofs; ++j) 
          {
            double val = bsetEn.evaluateSingle(j, volQuad, l, rhs_ );//*intel;
            singleRhs[j] += val;
          }
        }
        
        if(gradProblem_.hasSource())
        {
          // call source of gradient discrete model 
          gradCaller_.source(en, volQuad, l, gradSource_ );
        }

        if(gradProblem_.hasFlux())
        {
          // call anayltical flux of gradient discrete model 
          gradCaller_.analyticalFlux(en, volQuad, l, gradFMat_ );
        }

        if(problem_.hasFlux())
        {
          // call anayltical flux of discrete model 
          caller_.analyticalFlux(en, volQuad, l, fMat_ );
        }

        for(int k = 0; k < gradientNumDofs; ++k)
        {
          // eval tau_k 
          grdbsetEn.eval(k, volQuad, l, tau_[0] );
          
          if(gradProblem_.hasSource())
          {
            gradSourceNb_ = gradSource_; 
            // multiply tau with source
            for(int i=0; i<GradDimRange; ++i) gradSourceNb_[i] *= tau_[0][i];
            
            double val = grdbsetEn.evaluateSingle(k, volQuad, l, gradSourceNb_ ) * intel;

            // scalar product of basis functions is 1 (supposed to)
            massMatrixEn.add(k,k,val);
          }

          for (int j = 0; j < numDofs; ++j) 
          {
            // eval phi_j 
            bsetEn.eval(j, volQuad, l, phi_);
           
            helpmatr = gradFMat_; 
            for(int i = 0; i < dimDomain; ++i) 
            { 
              helpmatr[i][i] *= phi_[0];
            }

            {
              JacobianRangeType tmp;
              for(int i=0; i<dimGradRange; ++i) tmp[0][i] = tau_[0][i]*fMat_[0][i];
              
              // (tau_k)_l  * (grad phi_i)_l 
              double valC = bsetEn.evaluateGradientSingle(j, en, volQuad, l, tmp)*intel;       
              divMatrixEn.add( j , k , valC );
            }

            {
              // (phi_i)_l  * (div tau_k)_l
              double valB = -grdbsetEn.evaluateGradientSingle(k,en,volQuad,l,helpmatr)*intel;
              gradMatrixEn.add( k , j , valB );
            }
          }
        }
      } // end element integral 
  
      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
        // set default values 
        gradSource_   = 1.0;
        gradSourceNb_ = 1.0;

        int twistSelf = twistUtil_.twistInSelf(nit); 
        FaceQuadratureType faceQuadInner(nit, faceQuadOrd_, twistSelf, 
           FaceQuadratureType::INSIDE);
      
        // if neighbor exists 
        if (nit.neighbor()) 
        {
          EntityPointerType neighEp = nit.outside();
          EntityType&            nb = *neighEp;

          int twistNeighbor = twistUtil_.twistInNeighbor(nit);
          FaceQuadratureType faceQuadOuter(nit, faceQuadOrd_, twistNeighbor,
                                           FaceQuadratureType::OUTSIDE);
          
          caller_.setNeighbor(nb);
          gradCaller_.setNeighbor(nb);

          // create matrix handles for neighbor 
          MatrixAddHandleType stabMatrixNb(matrixHandler_.stabMatrix(),
                                         en, spc_, nb, spc_ ); 
          MatrixAddHandleType gradMatrixNb(matrixHandler_.gradMatrix(),
                                         en, gradientSpace_, nb, spc_ ); 
          MatrixAddHandleType divMatrixNb (matrixHandler_.divMatrix(),
                                         en , spc_, nb , gradientSpace_ ); 
          
          const BaseFunctionSetType& bsetNeigh = 
            spc_.getBaseFunctionSet(nb);
          const GradientBaseFunctionSetType& gradbsetNeigh = gradientSpace_.getBaseFunctionSet(nb);
         
          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l) 
          {
            DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
            double faceVol = unitNormal.two_norm();
            unitNormal *= 1.0/faceVol; 

            const double innerIntel = faceQuadInner.weight(l) * massVolElInv * faceVol ; 
            const double outerIntel = faceQuadOuter.weight(l) * massVolElInv * faceVol ; 

            if(gradProblem_.hasSource())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              gradCaller_.source(en, volQuad, l , gradSource_ );
              gradCaller_.source(nb, volQuad, l , gradSourceNb_ );
            }

            if(gradProblem_.hasFlux())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              gradCaller_.analyticalFlux(en, volQuad, l , gradFMat_ );
              gradCaller_.analyticalFlux(nb, volQuad, l , gradFMatNb_ );

              for(int i=0; i<dimGradRange; ++i) gradSource_[i] *= gradFMat_[0][i];
              for(int i=0; i<dimGradRange; ++i) gradSourceNb_[i] *= gradFMatNb_[0][i];
            }

            if(problem_.hasFlux())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              caller_.analyticalFlux(en, volQuad, l , fMat_   );
              caller_.analyticalFlux(nb, volQuad, l , fMatNb_ );

              for(int i=0; i<dimGradRange; ++i) gradSource_[i] *= fMat_[0][i];
              for(int i=0; i<dimGradRange; ++i) gradSourceNb_[i] *= fMatNb_[0][i];
            }

            ////////////////////////////////////////////////////////////
            //  
            //  FLUX evaluation 
            //  
            ////////////////////////////////////////////////////////////
            GradientRangeType sL,sR;
            RangeType uL,uR;

            // evaluate uFlux 
            gradCaller_.numericalFlux(nit,
                                      faceQuadInner,
                                      faceQuadOuter,
                                      l,
                                      sL,sR,
                                      uL,uR);

            RangeType enflux(uL);
            RangeType neighflux(uR);

            GradientRangeType sigmaEn,sigmaNb;

            // evaluate sigmaFlux 
            caller_.numericalFlux(nit,
                                  faceQuadInner,
                                  faceQuadOuter,
                                  l,
                                  sigmaEn,sigmaNb,
                                  uL,uR);

            double staben = uL; 
            double stabneigh = uR; 
            // to be revised 
            staben *= gradSource_[0]; 
            stabneigh *= gradSourceNb_[0];
            ////////////////////////////////////////////////////////////

            for(int i=0;i < gradientNumDofs;++i)
            {
              grdbsetEn.eval(i,faceQuadInner,l, tau_[0]);      
              gradbsetNeigh.eval(i,faceQuadOuter,l, tauneigh_[0]);
              
              RangeType valEn(0.0),valNeigh(0.0); 

              for(int j=0; j<numDofs; ++j)
              {
                //gradMAtrix
                bsetEn.   eval(j, faceQuadInner, l, phi_); 
                bsetNeigh.eval(j, faceQuadOuter, l, phiNeigh_ );
                
                {
                  double enVal= tau_[0] * unitNormal;
                  double neighVal= enVal;
                  
                  enVal    *=innerIntel;
                  enVal    *=enflux[0]*phi_[0];
                  
                  neighVal *=outerIntel;
                  neighVal *=neighflux[0]*phiNeigh_[0];
                 
                  // add value to matrix for en,nb
                  gradMatrixEn.add( i, j, enVal );

                  // add value to matrix for en,nb 
                  gradMatrixNb.add( i, j, neighVal );
                }

                {
                  double divmatValen,divmatValnb;
                  divmatValen = sigmaEn * tau_[0];
                  divmatValen *=phi_[0];
                  
                  divmatValen *=-innerIntel;
                   
                  divmatValnb = sigmaNb * tauneigh_[0];
                  divmatValnb *=phi_[0]; 
                  
                  divmatValnb *=-outerIntel; 
                
                  // add value to div matrix for en,en 
                  divMatrixEn.add( j, i , divmatValen );
                  
                  // add value to div matrix for en,nb 
                  divMatrixNb.add( j, i , divmatValnb );
                }

                if(i==0)
                {
                  for(int k=0;k<numDofs;++k)
                  {
                    double stabvalen = bsetEn.evaluateSingle(k, faceQuadInner, l,staben)    * innerIntel; 

                    double stabvalnb = bsetEn.evaluateSingle(k, faceQuadInner, l,stabneigh) * outerIntel;

                    // todo: make it right 
                    stabvalen *= phi_[0]; 

                    // add value to stqab matrix for en,en 
                    stabMatrixEn.add( k , j , stabvalen );

                    stabvalnb *= phiNeigh_[0];
                    // add value to stqab matrix for en,nb 
                    stabMatrixNb.add( k , j , stabvalnb );
                  }
                }
              }
            }
          }
        } // end if neighbor 

        // if intersection with boundary 
        if (nit.boundary()) 
        { 
          fMat_ = 1.0;
          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l) 
          {
            //const DomainType integrationNormal = nit.integrationOuterNormal(faceQuadInner.localPoint(l));
            DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
            double faceVol = unitNormal.two_norm();
            unitNormal *= 1.0/faceVol;
            
            const double bndFactor = faceQuadInner.weight(l) * massVolElInv;
            const double intel = bndFactor * faceVol;
            
            double t = 0.0;
            // get boundary value 
            RangeType boundaryValue(0.0);

            BoundaryIdentifierType bndType = problem_.boundaryValue(nit,t,
                faceQuadInner.localPoint(l),boundaryValue);

            if(gradProblem_.hasSource())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              gradCaller_.source(en, volQuad, l , gradSource_ );
            }

            if(gradProblem_.hasFlux())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              gradCaller_.analyticalFlux(en, volQuad, l , gradFMat_ );
              for(int i=0; i<dimGradRange; ++i) gradSource_[i] *= gradFMat_[0][i];
            }
            
            if(problem_.hasFlux())
            {
              // call anayltical flux 
              // todo: quadrature is not right here  
              caller_.analyticalFlux(en, volQuad, l , fMat_ );
              for(int i=0; i<dimGradRange; ++i) gradSource_[i] *= fMat_[0][i];
            }

            {
              RangeType fluxEn;
              GradientRangeType sigmaFluxEn,sigmaFluxFake; 
              caller_.boundaryFlux(nit, // intersection iterator 
                                   faceQuadInner,l, // quad and point number 
                                   sigmaFluxEn, fluxEn );
                
              for(int i=0; i<gradientNumDofs; ++i)
              { 
                grdbsetEn.eval(i,faceQuadInner,l, tau_[0]);

                // factor 2, becasue on boundary flux is identity, see
                // sigmaflux  
                double sigmaFlux = sigmaFluxEn * tau_[0];
                sigmaFlux *= 2.0 * intel;

                // dirichlet boundary value for sigma 
                if(bndType.isDirichletNonZero())
                {
                  // u^ on \gamma_D
                  RangeType bndVal(boundaryValue);
                  bndVal *= tau_[0] * unitNormal; // sigmaFlux; 
                  bndVal *= intel;
                  //bndVal *= gradFMat_[0][0];
                  gradRhs[i] += bndVal[0];
                }
                              
                for (int j = 0; j < numDofs; ++j) 
                {
                  bsetEn.eval(j, faceQuadInner, l, phi_);
                 
                  // u+ in u^
                  double Val = sigmaFlux * phi_[0];
           
                  if(bndType.isDirichletType())
                  {
                    Val*= fMat_[0][0];
                    // add value for dirichlet bnd 
                    divMatrixEn.add( j , i , -Val );
                  }

                  if(bndType.isNeumannType())
                  {
                    Val *= gradFMat_[0][0];
                    gradMatrixEn.add( i , j , Val );
                  }
                  
                  if(bndType.isDirichletType())
                  {
                    if(i==0)
                    {
                      // dirichlet boundary values for u 
                      if(bndType.isDirichletNonZero())
                      {
                        RangeType bndVal (boundaryValue);
                          
                        // only valid for dim range = 1
                        double rhsVal1 = bndVal[0] * phi_[0];
                        
                        rhsVal1 *= bndFactor;
                        rhsVal1 *= gradSource_[0];
                        singleRhs[j] += rhsVal1;
                      }
                      
                      for(int k=0; k<numDofs; ++k)
                      {
                        double stabvalen=0.0;
                        
                        // note that only gLeft is used here 
                        RangeType fluxEnTmp =  fluxEn[0] * phi_[0];  
                        
                        stabvalen = bsetEn.evaluateSingle(k, faceQuadInner, l,fluxEnTmp) * intel;
                        
                        stabvalen *= gradSource_[0];
                        //stabvalen *= fMat_[0][0];
                        // add value to stab matrix for en,en 
                        stabMatrixEn.add( k , j , stabvalen );
                      }
                    }
                  }

                  if(bndType.isNeumannNonZero())
                  {
                    // Neumann boundary
                    if(i==0)
                    {
                      RangeType bndVal (boundaryValue);
                      // dirichlet boundary values for u 
                      bndVal *= phi_[0];
                      bndVal *= bndFactor;
                      singleRhs[j] += bndVal[0];
                    }
                  }
                } 
              }
            }
          }
        } // end if boundary

      } // end intersection iterator 

    } // end apply local 

    void update(const ArgumentType& arg, DestinationType& dest) const
    {
      prepare(arg, dest);

      typedef typename DiscreteFunctionSpaceType :: IteratorType  IteratorType ;
      IteratorType endit = spc_.end();
      for (IteratorType it = spc_.begin(); it != endit; ++it)
      {
        updateLocal(*it);
      }

      finalize(arg, dest);
    }

    void updateMatrix( const ArgumentType & arg,
                       DestinationType & rhs )
    {
      if(matrixHandler_.hasMassMatrix())
      {
        matrixHandler_.clearMass();

        this->update(arg,rhs);
        matrixAssembled_ = true;

        createPreconditionMatrix();
      }
    }

    // calculate pre-condition matrix 
    void createPreconditionMatrix()
    {
      if(matrixHandler_.hasPcMatrix())
      {
        const int singleSize = spc_.size();

        assert( diag_ );
        DestinationType & diag = *diag_; 
        
        if(gradProblem_.hasSource())
        {
          matrixHandler_.divMatrix().getDiag( matrixHandler_.massMatrix(), matrixHandler_.gradMatrix() , diag );
        }
        else 
        {
          matrixHandler_.divMatrix().getDiag( matrixHandler_.gradMatrix() , diag );
        }

        matrixHandler_.stabMatrix().addDiag( diag );

        double * diagPtr = diag.leakPointer();
        for(register int i=0; i<singleSize; ++i) 
        {
          double val = diagPtr[i]; 
          // when using parallel Version , we could have zero on diagonal
          // for ghost elements 
          assert( (spc_.grid().comm().size() > 1) ? 1 : (std::abs( val ) > 0.0 ) );
          if( std::abs( val ) > 0.0 )
          {
            val = 1.0/val; 
            matrixHandler_.pcMatrix().add(i,i,val);
          }
        }
      }
    }

    void updateLocal(EntityType& en) const
    {
      // local function for gradient right hand side 
      typedef typename GradDestinationType :: LocalFunctionType GradLFType; 
      GradLFType gradRhs = gradRhs_.localFunction(en); //rhs
      
      //- typedefs
      typedef typename DiscreteFunctionSpaceType::IndexSetType IndexSetType;
      typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

      MatrixAddHandleType massMatrixEn (matrixHandler_.massMatrix(),
                                     en, gradientSpace_, en, gradientSpace_ ); 
      
      typedef typename DiscreteGradientSpaceType::IndexSetType GradientIndexSetType;
      typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;

      //- statements
      caller_.setEntity(en);
      gradCaller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      const VolumeQuadratureType massQuad(en, 0 );   
      const GeometryType & geo = en.geometry();
      const double massVolElInv = massVolumeInv(geo,massQuad);
      
      const GradientBaseFunctionSetType& grdbsetEn = gradientSpace_.getBaseFunctionSet(en);
      const int gradientNumDofs = grdbsetEn.numBaseFunctions();
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////

      //set default values 
      gradSource_   = 1.0;
      gradSourceNb_ = 1.0;
      // set to id matrix 
      gradFMat_   = one_;
      gradFMatNb_ = one_;

      fMat_       = 1.0;
      fMatNb_     = 1.0;
      
      const int quadNop = volQuad.nop();
      GradJacobianRangeType helpmatr(one_);

      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc factor for bas functions 
        double intel = volQuad.weight(l)*
            geo.integrationElement(volQuad.point(l))*massVolElInv;
        
        if(gradProblem_.hasSource())
        {
          // call source of gradient discrete model 
          gradCaller_.source(en, volQuad, l, gradSource_ );
        }

        for(int k = 0; k < gradientNumDofs; ++k)
        {
          // eval tau_k 
          grdbsetEn.eval(k, volQuad, l, tau_[0] );
          
          if(gradProblem_.hasSource())
          {
            gradSourceNb_ = gradSource_; 
            // multiply tau with source
            for(int i=0; i<GradDimRange; ++i) gradSourceNb_[i] *= tau_[0][i];
            
            double val = grdbsetEn.evaluateSingle(k, volQuad, l, gradSourceNb_ ) * intel;

            // scalar product of basis functions is 1 (supposed to)
            massMatrixEn.add(k,k,val);
          }
        }
      } // end element integral 
    }
  
  private:
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    //! empty constructor not defined 
    LocalDGElliptOperator();
    //! copy constructor not defined 
    LocalDGElliptOperator(const LocalDGElliptOperator&);

  private:
    double massVolumeInv(const GeometryType& geo, const VolumeQuadratureType & quad ) const
                         
    {
      double result = 0.0;
      double massVolInv = 0.0;
      const int quadNop = quad.nop();
      for (int qp = 0; qp < quadNop; ++qp) 
      {
        massVolInv += quad.weight(qp);//volumen referenzelement
        result += 
          quad.weight(qp) * geo.integrationElement(quad.point(qp));
      }
      massVolInv /= result;
      return massVolInv;
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    GradientPassType & gradPass_; 
    mutable GradientModelCallerType & gradCaller_; 
    GradientDiscreteModelType& gradProblem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const DiscreteGradientSpaceType & gradientSpace_;
    mutable SingleCommunicationManagerType singleCommunicate_;
    mutable GradCommunicationManagerType gradCommunicate_;
    
    mutable GradDestinationType gradTmp_; 
    mutable GradDestinationType gradRhs_; 
    mutable GradDestinationType * massTmp_; 
    mutable DestinationType multTmp_;
    mutable DestinationType * diag_;
    mutable double dtMin_;
  
    //! Some helper variables
    mutable JacobianRangeType fMat_;
    mutable JacobianRangeType fMatNb_;
    mutable GradientRangeType gradSource_;
    mutable GradientRangeType gradSourceNb_;
    mutable GradientJacobianRangeType one_;
    mutable GradientJacobianRangeType gradFMat_;
    mutable GradientJacobianRangeType gradFMatNb_;
    mutable RangeType valEn_;
    mutable RangeType valNeigh_;
    mutable RangeType baseEn_;
    mutable RangeType baseNeigh_;
    mutable RangeType phi_;
    mutable RangeType source_;
    mutable RangeType sourceNb_;
    mutable RangeType rhs_;

    //mutable GradientRangeType tau_;
    mutable JacobianRangeType tau_;
    mutable JacobianRangeType tauneigh_;
    mutable RangeType phiNeigh_;
    mutable RangeType gradEval_;
    mutable DomainType grads_;

    TimeProvider* time_;

    TwistUtility<GridType> twistUtil_;

    int volumeQuadOrd_,faceQuadOrd_;

    mutable MatrixHandlerType matrixHandler_;

    mutable bool matrixAssembled_;
    mutable bool preCond_;
  };
  
  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGElliptPass :
    public LocalPass<DiscreteModelImp, typename PreviousPassImp:: PreviousPassType> 
  {
  public:
    typedef typename PreviousPassImp::PreviousPassType  PreviousPassType;
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassType> BaseType;

    typedef LocalDGElliptPass<DiscreteModelImp,PreviousPassImp> ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp GradFePassImp;

    typedef typename GradFePassImp :: DestinationType GradDestinationType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType;
    typedef typename EntityType :: EntityPointer EntityPointerType;
    
    typedef typename BaseType::ArgumentType ArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;

    // Types extracted from the underlying grids
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;
   
    // Range of the destination
    enum { dimR = DiscreteFunctionSpaceType::DimRange };
    enum { dimD = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange  = DiscreteFunctionSpaceType :: DimRange }; 
    enum { dimDomain = DiscreteFunctionSpaceType :: DimDomain }; 
                    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;

    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };

    // define previous pass of grad pass as previous pass of hole ellipt
    // pass
    typedef typename GradFePassImp :: PreviousPassType ElliptPrevPassType;
    // define ellipt operator 
    typedef LocalDGElliptOperator<DiscreteModelImp,GradFePassImp,ElliptPrevPassType> FEOperatorType;

    DiscreteModelType& problem_; 
    DiscreteFunctionSpaceType& spc_;
    
    //ElliptPrevPassType feStartPass_;
    mutable FEOperatorType op_;

    // define type of inverse operator 
    typedef typename DiscreteModelType :: Traits :: Traits :: template 
      InverseOperator<DestinationType,FEOperatorType>:: InverseOperatorType InverseOperatorType;

    const double eps_;
    const int maxIterFactor_; 
    mutable int maxIter_;

    InverseOperatorType invOp_; 

    mutable DestinationType rhs_;

    mutable int sequence_;
    const bool verbose_;
      
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    LocalDGElliptPass(DiscreteModelType& problem, 
                PreviousPassImp & pass, 
                DiscreteFunctionSpaceType& spc,
                double eps = 1e-10 , int maxIterFactor = 3 , bool verbose = false )
      : BaseType(pass.previousPass(),spc)
      , problem_(problem)
      , spc_(spc) 
      , op_(problem,pass,pass.previousPass(),spc)
      , eps_(eps)
      , maxIterFactor_(maxIterFactor) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps,eps,maxIter_,verbose)
      , rhs_("FEPass::RHS",spc)
      , sequence_(-1)
      , verbose_(verbose)
    {}

    void applyLocal(EntityType& en) const
    {
    }

    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // if grid has changed, recalc matrices  
      if(sequence_ != spc_.sequence())
      {
        // only clear destination if matrix has really changed 
        // otherwise keep old value as initial value 
        dest.clear();
        op_.prepare(arg,rhs_);
        op_.buildMatrix(arg,  rhs_ );
        sequence_ = spc_.sequence();
      }
      else 
      {
        op_.updateMatrix( arg, rhs_ );
      }
    }

    //! Some timestep size management.
    virtual void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    virtual void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      // prepare operator 
      prepare(arg,dest);

      // calculate new maxIter  
      maxIter_ = maxIterFactor_ * spc_.size();

      // set parameter which might have changed
      //invOp_.setParameters(eps_,maxIter_,verbose_);

      // solve the system 
      invOp_(rhs_,dest);
    } 

    template <class FuncType, class GradType>
    void evalGradient(const FuncType & u, GradType & grad ) const
    {
      op_.evalGradient(u,grad);
    }
  };

  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class PreviousPassImp>
  class LocalDGElliptGradPass :
    public LocalPass<DiscreteModelImp, PreviousPassImp> 
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    typedef LocalDGElliptGradPass<DiscreteModelImp,PreviousPassImp> ThisType;

    //! Repetition of template arguments
    typedef DiscreteModelImp DiscreteModelType;
    //! Repetition of template arguments
    typedef PreviousPassImp PreviousPassType;

    // Types from the base class
    typedef typename BaseType::Entity EntityType; 
    typedef typename EntityType::EntityPointer EntityPointerType;
    typedef typename BaseType::ArgumentType ArgumentType;
    typedef typename BaseType::GlobalArgumentType GlobalArgumentType;

    // Types from the traits
    typedef typename DiscreteModelType::Traits::DestinationType DestinationType;
    typedef typename DiscreteModelType::Traits::VolumeQuadratureType VolumeQuadratureType;
    typedef typename DiscreteModelType::Traits::FaceQuadratureType FaceQuadratureType;
    typedef typename DiscreteModelType::Traits::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
        
    
    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridType GridType;
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    
    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry Geometry;


    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    
    //my Typedefs
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd =DiscreteFunctionSpaceType::polOrd};

    typedef typename DiscreteFunctionSpaceType::DomainFieldType DomainFieldType;
    typedef typename DiscreteFunctionSpaceType::RangeFieldType RangeFieldType;

  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param quadOrd0 defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param quadOrd1 defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGElliptGradPass(DiscreteModelType& problem, 
                    PreviousPassType& pass, 
                    DiscreteFunctionSpaceType& spc) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      spc_(spc),
      prevPass_(pass)
    {
    }

    //! don't allocate memory here 
    virtual void allocateLocalMemory() {}
   
    //! return reference to caller 
    DiscreteModelCallerType & caller () { return caller_; }

    //! return previous pass of this pass 
    PreviousPassType & previousPass() { return this->previousPass_; }

    //! return problem for real fe pass 
    DiscreteModelType & problem () { return problem_; }

    const DiscreteFunctionSpaceType & space () const { return spc_; }
   
    //! Destructor
    virtual ~LocalDGElliptGradPass() {
    }

    void applyLocal(EntityType& en) const
    {
    }

    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      // normal call procedure 
      prevPass_.pass(arg);

      // now get gradient from previous pass 
      prevPass_.evalGradient(prevPass_.destination(),dest);
    }

    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! Some timestep size management.
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    DiscreteFunctionSpaceType& spc_;
    mutable PreviousPassImp & prevPass_;
  };

} // end namespace Dune

#endif
