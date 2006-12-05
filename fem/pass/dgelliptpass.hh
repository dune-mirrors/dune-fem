#ifndef DUNE_ELLIPTPASS_HH
#define DUNE_ELLIPTPASS_HH

//- system includes 
#include <set>

//- Dune includes 
#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/utility/twistutility.hh>

//- local includes 
#include <dune/fem/pass/pass.hh>
#include <dune/fem/pass/discretemodel.hh>
#include <dune/fem/pass/modelcaller.hh>

#include <dune/fem/misc/timeutility.hh>
#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/discretefunction/common/dfcommunication.hh>
#include <dune/fem/space/common/communicationmanager.hh>

namespace Dune {
/*! @defgroup PassEllipt Local Discontinous Galerkin for second order elliptic equations
 *  @ingroup Pass
 * Description: Solver for equations of the form
** \f{eqnarray*}
**   div(A(x)\nabla u) + &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/

  //! Concrete implementation of Pass for elliptic equations using LDG
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

    //! return reference to space 
    const DiscreteFunctionSpaceType & space () const { return spc_; }
   
    //! Destructor
    virtual ~LocalDGElliptGradientPass() {
    }

    //! do nothing here 
    void applyLocal(EntityType& en) const
    {
    }

    //! do nothing here  
    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      abort();
    }

    //! do nothing here 
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
    typedef typename GridPartType :: GridType :: Traits :: LocalIdSet LocalIdSetType; 
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
    typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;

    typedef typename DiscreteGradientSpaceType::RangeType GradientRangeType;
    typedef typename DiscreteGradientSpaceType::JacobianRangeType GradJacobianRangeType;
    typedef GradJacobianRangeType GradientJacobianRangeType;
    enum { GradDimRange = GradientRangeType :: dimension };
    
    typedef typename DiscreteGradientSpaceType::BaseFunctionSetType GradientBaseFunctionSetType;
    
    typedef typename DiscreteModelType :: Traits :: Traits ::template
      MatrixHandler<DiscreteFunctionSpaceType,DiscreteGradientSpaceType> ::
      MatrixHandlerType MatrixHandlerType; 

    typedef typename MatrixHandlerType::MatrixAddHandleType MatrixAddHandleType;
    typedef typename MatrixHandlerType::MatrixType MatrixType;
    typedef typename MatrixHandlerType::PreconditionMatrixType PreconditionMatrixType;
    
    typedef typename DiscreteModelType :: BoundaryIdentifierType BoundaryIdentifierType;    

    typedef CommunicationManager<DiscreteFunctionSpaceType> SingleCommunicationManagerType; 
    typedef CommunicationManager<DiscreteGradientSpaceType> GradCommunicationManagerType; 

    typedef typename LocalIdSetType :: IdType LocalIdType;
    typedef std::set< LocalIdType > EntityMarkerType; 
    //typedef std::map< LocalIdType , bool > EntityMarkerType; 
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGElliptOperator(DiscreteModelType& problem, 
                GradientPassType & gradPass,
                PreviousPassType& pass, 
                DiscreteFunctionSpaceType& spc,
                bool verbose = false, 
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
      localIdSet_(gridPart_.grid().localIdSet()),
      gradientSpace_(gradPass_.space()),
      singleCommunicate_(spc_),
      gradCommunicate_(gradientSpace_),
      gradTmp_("FEPass::gradTmp",gradientSpace_),
      gradRhs_("FEPass::gradRhs",gradientSpace_),
      massTmp_(0),
      //multTmp_("FEPass::multTmp",spc_),
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
      elemOrder_(std::max(spc_.order(),gradientSpace_.order())),
      faceOrder_(std::max(spc_.order(),gradientSpace_.order())+1),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (2*elemOrder_) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (2*faceOrder_) : faceQuadOrd ),
      matrixHandler_(spc_,gradientSpace_,gradProblem_.hasSource(),problem_.preconditioning()),
      entityMarker_(),
      matrixAssembled_(false),
      verbose_(verbose)
    {
      assert( matrixHandler_.hasMassMatrix() == gradProblem_.hasSource() );
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );

      // need for multTmpPointer 
      assert( spc_.size() <= gradientSpace_.size() );

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
        std::cerr << "Source for DGElliptPass not supported yet. \n";
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

    //! setup matrix 
    void buildMatrix(const ArgumentType & arg, DestinationType & rhs)
    {
      // resize matrices 
      matrixHandler_.resize(verbose_);

      // clear matrices 
      matrixHandler_.clear();
     
      gradRhs_.clear();
      rhs.clear();
      
      dest_ = &rhs;

      // build matrix and rhs 
      this->compute( arg, rhs );
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

      // adjust rhs 
      double * multTmpPointer = gradTmp_.leakPointer();
      matrixHandler_.divMatrix().multOEM(rhsPtr, multTmpPointer);
      
      double * singleRhsPtr = rhs.leakPointer();
      const int singleSize = spc_.size();
      for(register int i=0; i<singleSize; ++i) 
      {
        singleRhsPtr[i] -= multTmpPointer[i];      
      }
    }

    //! rebuild matrix after adaptation 
    void reBuildMatrix(const ArgumentType & arg, DestinationType & rhs,
        bool check )
    {
      if( check && entityMarker_.empty() ) 
      {
        buildMatrix( arg, rhs );
        return ;
      } 
      
      // set dest to rhs 
      dest_ = &rhs;

      // resize matrices 
      matrixHandler_.resize(verbose_);

      // prepare operator 
      prepare (arg, rhs);

      // marking should be done in restrictLocal and prolongLocal
      // find better way to also makr all neighbors of marked entities 
      typedef typename DiscreteFunctionSpaceType::IteratorType IteratorType;
      {
        IteratorType end = spc_.end();
        for(IteratorType it = spc_.begin(); it != end; ++it)
        {
          EntityType & en = *it; 
          if(en.wasRefined() || ( entityMarker_.find( localIdSet_.id( en ) ) != entityMarker_.end() ) )
          {
            entityMarker_.insert( localIdSet_.id( en ) );
            IntersectionIteratorType endnit = gridPart_.iend(en); 
            for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
            { 
              if(nit.neighbor())
              {
                EntityPointerType ep = nit.outside();
                EntityType & nb = *ep;
                entityMarker_.insert( localIdSet_.id( nb ) );
              }
            }
          }
        }
      }

      // vector to new all new dofs 
      std::vector<int> newDofs;

      // if only one geometry type, we can reserve memory at once 
      if( ! spc_.multipleGeometryTypes() )
      {
        IteratorType it = spc_.begin();
        if( it != spc_.end() )
        {
          // get base function set of single space 
          typedef typename DiscreteFunctionSpaceType::BaseFunctionSetType BaseFunctionSetType;
          const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(*it);
          const int numDofs = bsetEn.numBaseFunctions();
          newDofs.reserve( numDofs * entityMarker_.size() );
        }
      }

      // reset vector 
      newDofs.resize( 0 );

      // rebuild matrix for all marked entities 
      {
        IteratorType end = spc_.end();
        typedef typename EntityMarkerType :: iterator iterator ;
        iterator markerEnd = entityMarker_.end();
        for(IteratorType it = spc_.begin(); it != end; ++it)
        {
          EntityType & en = *it; 
          // if entity was marked for rebuild 
          if( entityMarker_.find( localIdSet_.id(en ) ) != markerEnd )
          {
            // clear all entries belonging to this entity 
            // and create map for all dofs that belong to spc_ for this entity
            clearLocal(en, newDofs );

            // build local matrices 
            applyLocal(en);
          }
        }
      }

      // unset pointer to discrete functions 
      finalize( arg, rhs );

      // matrix is build up now 
      matrixAssembled_ = true;

      // create pre-condition matrix if activated 
      createPreconditionMatrix();

      // adjust right hand side 
      double * rhsPtr = gradRhs_.leakPointer();
      if(gradProblem_.hasSource())
      {
        double * massTmpPointer = massTmp_->leakPointer();
        matrixHandler_.massMatrix().multOEM(rhsPtr,massTmpPointer);
        rhsPtr = massTmpPointer;
      } 

      // adjust rhs 
      double * multTmpPointer = gradTmp_.leakPointer();
      matrixHandler_.divMatrix().multOEM(rhsPtr,multTmpPointer );

      double * singleRhsPtr = rhs.leakPointer();

      // adjust only new vlaues of rhs 
      // ( diag contians 0.0 where old values are and 1.0 where new )
      const int size = newDofs.size();
      for(int i=0; i<size; ++i) 
      {
        singleRhsPtr[ newDofs[i] ] -= multTmpPointer[ newDofs[i] ];      
      }

      // empty all marked entities 
      entityMarker_.clear();
    }

    //! calculates M^-1 B * u, is u is the solution then outcome is grad u 
    template <class FuncType, class GradType> 
    void evalGradient(const FuncType & u, GradType & grad) const
    {
      grad.clear();
      // if source then apply also mass matrix 
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
        // only apply grad matrix here 
        matrixHandler_.gradMatrix().multOEM(u.leakPointer(),grad.leakPointer());
        grad += gradRhs_;
      }
    }

    //! do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer());
    }
    
    //! do matrix vector multiplication, used by OEM-Solver and DuneODE Solvers  
    void multOEM(const double * arg, double * dest) const
    {
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
        gradCommunicate_.exchange( *massTmp_ );
      }
      else 
      {
        // otherwise communicate grad 
        gradCommunicate_.exchange( gradTmp_ );
      }

      // calc dest = divMatrix * gradTmp 
      matrixHandler_.divMatrix().multOEM(gradTmpPointer, dest );
      // calc dest += stabMatrix * arg 
      matrixHandler_.stabMatrix().multOEMAdd(arg,dest);
    }

  private:   
    //! create discrete function from double * and communicate data 
    //! assumes that x has size of single space 
    void communicate(const double * x) const
    {
      DestinationType dest("DGEllipt::communicate_tmp",spc_,x);
      singleCommunicate_.exchange( dest );
    }
      
  public:
    //! return refernence to system matrix, used by OEM-Solver
    const ThisType & systemMatrix () const { return *this; }
    
    //! return reference to preconditioning matrix, used by OEM-Solver
    const PreconditionMatrixType & preconditionMatrix () const { return matrixHandler_.pcMatrix(); }

    //! returns true if preconditioning matrix has been build 
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
    //! return reference to function space 
    const DiscreteFunctionSpaceType & space() const 
    {
      return spc_;
    }

  public:
    //! set all data belonging to this entity to zero 
    void restrictLocal(const EntityType& father, const EntityType & son, bool firstCall ) const
    {
      assert( gridPart_.indexSet().adaptive() );
      if( firstCall ) 
      {
        // insert new entity 
        entityMarker_.insert( localIdSet_.id( father ) );
      }
    }

    //! set all data belonging to this entity to zero 
    void prolongLocal(const EntityType& father, const EntityType & son, bool ) const
    {
      assert( gridPart_.indexSet().adaptive() );
      // insert new entity  
      entityMarker_.insert( localIdSet_.id( son ) );
    }

    //! set all data belonging to this entity to zero 
    void clearLocal(EntityType& en, std::vector<int> & newDofs ) const
    {
      {
        // local function for right hand side 
        typedef typename DestinationType :: LocalFunctionType SingleLFType; 
        SingleLFType singleRhs = dest_->localFunction(en); //rhs

        const int numDofs = singleRhs.numDofs();
        if( spc_.multipleGeometryTypes() )
        {
          newDofs.reserve( newDofs.size() + numDofs );
        }

        assert( numDofs > 0 );
        for(int i=0; i<numDofs; ++i) 
        {
          singleRhs[i] = 0.0;
          // remember dof num of new dofs 
          newDofs.push_back( spc_.mapToGlobal( en, i ) );
        }
      }
      
      {
        // local function for gradient right hand side 
        typedef typename GradDestinationType :: LocalFunctionType GradLFType; 
        GradLFType gradRhs = gradRhs_.localFunction(en); //rhs

        const int numDofs = gradRhs.numDofs();
        for(int i=0; i<numDofs; ++i) gradRhs[i] = 0.0;
      }
      
      {
        MatrixAddHandleType stabMatrixEn(matrixHandler_.stabMatrix(),
                                       en, spc_, en, spc_ ); 
        stabMatrixEn.clear();
      }
      
      MatrixAddHandleType gradMatrixEn(matrixHandler_.gradMatrix(),
                                     en, gradientSpace_, en, spc_ ); 
      gradMatrixEn.clear();
      MatrixAddHandleType divMatrixEn (matrixHandler_.divMatrix(),
                                     en, spc_, en, gradientSpace_ ); 
      divMatrixEn.clear();
      
      if(gradProblem_.hasSource())
      {
        MatrixAddHandleType massMatrixEn (matrixHandler_.massMatrix(),
                                     en, gradientSpace_, en, gradientSpace_ ); 
        massMatrixEn.clear();
      }
    }
    
    //! apply operator on entity 
    void applyLocal(EntityType& en) const
    {
      // local function for right hand side 
      typedef typename DestinationType :: LocalFunctionType SingleLFType; 
      SingleLFType singleRhs = dest_->localFunction(en); //rhs
      
      // local function for gradient right hand side 
      typedef typename GradDestinationType :: LocalFunctionType GradLFType; 
      GradLFType gradRhs = gradRhs_.localFunction(en); //rhs
      
      MatrixAddHandleType stabMatrixEn(matrixHandler_.stabMatrix(),
                                     en, spc_, en, spc_ ); 
      MatrixAddHandleType gradMatrixEn(matrixHandler_.gradMatrix(),
                                     en, gradientSpace_, en, spc_ ); 
      MatrixAddHandleType divMatrixEn (matrixHandler_.divMatrix(),
                                     en, spc_, en, gradientSpace_ ); 
      
      MatrixAddHandleType massMatrixEn (matrixHandler_.massMatrix(),
                                     en, gradientSpace_, en, gradientSpace_ ); 
      
      // make entities known in callers
      caller_.setEntity(en);
      gradCaller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      const GeometryType & geo = en.geometry();
      const double massVolElInv = massVolumeInv(geo);
      
      // get base function set of single space 
      const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      
      // get base function set of gradient space 
      const GradientBaseFunctionSetType& grdbsetEn = gradientSpace_.baseFunctionSet(en);
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
        const double intel = volQuad.weight(l)*
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

        // if neighbor exists 
        if (nit.neighbor()) 
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
            applyLocalNeighbor(nit,en,massVolElInv,volQuad,
                  faceQuadInner,faceQuadOuter, 
                  bsetEn,grdbsetEn,
                  stabMatrixEn,gradMatrixEn,divMatrixEn);
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
            applyLocalNeighbor(nit,en,massVolElInv,volQuad,
                  nonConformingFaceQuadInner,
                  nonConformingFaceQuadOuter, 
                  bsetEn,grdbsetEn,
                  stabMatrixEn,gradMatrixEn,divMatrixEn);
          }
        } // end if neighbor 

        // if intersection with boundary 
        if (nit.boundary()) 
        { 
          fMat_ = 1.0;
          
          FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                           FaceQuadratureType::INSIDE);

          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l) 
          {
            //const DomainType integrationNormal = nit.integrationOuterNormal(faceQuadInner.localPoint(l));
            DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
            const double faceVol = unitNormal.two_norm();
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

            {
              RangeType sigmaFluxEn; 
              
              caller_.boundaryFlux(nit, // intersection iterator 
                                   faceQuadInner,l, // quad and point number 
                                   sigmaFluxEn);
                
              for(int i=0; i<gradientNumDofs; ++i)
              { 
                grdbsetEn.eval(i,faceQuadInner,l, tau_[0]);

                // factor 2, becasue on boundary flux is identity, see
                // sigmaflux  
                double sigmaFlux = unitNormal * tau_[0];
                sigmaFlux *= 2.0 * intel * sigmaFluxEn;

                // dirichlet boundary value for sigma 
                if(bndType.isDirichletNonZero())
                {
                  // u^ on \gamma_D
                  RangeType bndVal(boundaryValue);
                  bndVal *= tau_[0] * unitNormal; // sigmaFlux; 
                  bndVal *= intel;
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
                        RangeType fluxEnTmp = phi_[0];  
                        
                        stabvalen = bsetEn.evaluateSingle(k, faceQuadInner, l,fluxEnTmp) * bndFactor;
                        
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

      // resort corresponding matrix rows for ascending numbering 
      stabMatrixEn.resort(); 
      gradMatrixEn.resort();
      divMatrixEn.resort();
     
      if( gradProblem_.hasSource())
      {
        massMatrixEn.resort();
      }
    } // end apply local 


    template <class QuadratureImp> 
    void applyLocalNeighbor(IntersectionIteratorType & nit, 
                            EntityType & en, const double massVolElInv,
                            VolumeQuadratureType & volQuad,
                            const QuadratureImp & faceQuadInner, 
                            const QuadratureImp & faceQuadOuter, 
                            const BaseFunctionSetType & bsetEn, 
                            const GradientBaseFunctionSetType & grdbsetEn, 
                            MatrixAddHandleType & stabMatrixEn, 
                            MatrixAddHandleType & gradMatrixEn, 
                            MatrixAddHandleType & divMatrixEn) const
    {
      const int numDofs = bsetEn.numBaseFunctions();
      const int gradientNumDofs = grdbsetEn.numBaseFunctions();

      EntityPointerType neighEp = nit.outside();
      EntityType&            nb = *neighEp;

      caller_.setNeighbor(nb);
      gradCaller_.setNeighbor(nb);

      // create matrix handles for neighbor 
      MatrixAddHandleType stabMatrixNb(matrixHandler_.stabMatrix(),
                                     en, spc_, nb, spc_ ); 
      MatrixAddHandleType gradMatrixNb(matrixHandler_.gradMatrix(),
                                     en, gradientSpace_, nb, spc_ ); 
      MatrixAddHandleType divMatrixNb (matrixHandler_.divMatrix(),
                                     en , spc_, nb , gradientSpace_ ); 
      
      // get base function set 
      const BaseFunctionSetType& bsetNeigh = spc_.baseFunctionSet(nb);
      const GradientBaseFunctionSetType& gradbsetNeigh 
                              = gradientSpace_.baseFunctionSet(nb);
     
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
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

        ////////////////////////////////////////////////////////////
        //  
        //  FLUX evaluation 
        //  
        ////////////////////////////////////////////////////////////
        GradientRangeType uFluxSigmaLeft,uFluxSigmaRight; // sigma parts of uFlux 
        GradientRangeType uFluxULeft,uFluxURight; // u parts of uFlux 

        // evaluate uFlux, return type is vector  
        gradCaller_.numericalFlux(nit,
                                  faceQuadInner,
                                  faceQuadOuter,
                                  l,
                                  uFluxSigmaLeft,uFluxSigmaRight,
                                  uFluxULeft,uFluxURight);

        RangeType sigmaFluxSigmaLeft,sigmaFluxSigmaRight; // sigma  parts of sigmaFlux 
        RangeType sigmaFluxULeft,sigmaFluxURight; // u parts of sigmaFlux 

        // evaluate sigmaFlux, return type is scalar 
        caller_.numericalFlux(nit,
                              faceQuadInner,
                              faceQuadOuter,
                              l,
                              sigmaFluxSigmaLeft,sigmaFluxSigmaRight,
                              sigmaFluxULeft,sigmaFluxURight);

        // to be revised 
        sigmaFluxULeft  *= gradSource_[0]; 
        sigmaFluxURight *= gradSourceNb_[0];
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
              // eval tau * (un)  (scalar procduct) 
              double enVal    = tau_[0] * uFluxULeft;
              double neighVal = tau_[0] * uFluxURight; 
              
              enVal    *= innerIntel;
              enVal    *= phi_[0];
              
              neighVal *= outerIntel;
              neighVal *= phiNeigh_[0];
             
              // add value to matrix for en,nb
              gradMatrixEn.add( i, j, enVal );

              // add value to matrix for en,nb 
              gradMatrixNb.add( i, j, neighVal );
            }

            {
              // value en 
              double divmatValen = unitNormal * tau_[0];
              divmatValen *= sigmaFluxSigmaLeft;
              divmatValen *= phi_[0];
              divmatValen *=-innerIntel;
               
              // value neighbor 
              double divmatValnb = unitNormal * tauneigh_[0];
              divmatValnb *= sigmaFluxSigmaRight;
              divmatValnb *= phi_[0]; 
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
                double stabvalen = bsetEn.evaluateSingle(k, faceQuadInner, l, sigmaFluxULeft) 
                                   * innerIntel; 

                double stabvalnb = bsetEn.evaluateSingle(k, faceQuadInner, l, sigmaFluxURight) 
                                   * outerIntel;

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
    }
    //! only calculate mass matrix new
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

    //! only calculate mass matrix new on entity 
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
        matrixHandler_.clearPcMatrix();

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

      const GeometryType & geo = en.geometry();
      const double massVolElInv = massVolumeInv(geo);
      
      const GradientBaseFunctionSetType& grdbsetEn = gradientSpace_.baseFunctionSet(en);
      const int gradientNumDofs = grdbsetEn.numBaseFunctions();
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////

      //set default values 
      GradientRangeType gradSource(1.0);
      GradientRangeType gradSourceNb(1.0);
      
      // cache number of integration points 
      const int quadNop = volQuad.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc factor for bas functions 
        const double intel = volQuad.weight(l)*
            geo.integrationElement(volQuad.point(l))*massVolElInv;
        
        if(gradProblem_.hasSource())
        {
          // call source of gradient discrete model 
          gradCaller_.source(en, volQuad, l, gradSource );
        }

        for(int k = 0; k < gradientNumDofs; ++k)
        {
          // eval tau_k 
          grdbsetEn.eval(k, volQuad, l, tau_[0] );
          
          if(gradProblem_.hasSource())
          {
            gradSourceNb = gradSource; 
            // multiply tau with source
            for(int i=0; i<GradDimRange; ++i) gradSourceNb[i] *= tau_[0][i];
            
            double val = grdbsetEn.evaluateSingle(k, volQuad, l, gradSourceNb ) * intel;

            // scalar product of basis functions is 1 (supposed to)
            massMatrixEn.add(k,k,val);
          }
        }
      } // end element integral 
    }

    GradDestinationType & tmpMemory() { return gradTmp_; }
  
  private:
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    //! empty constructor not defined 
    LocalDGElliptOperator();
    //! copy constructor not defined 
    LocalDGElliptOperator(const LocalDGElliptOperator&);

  private:
    double massVolumeInv(const GeometryType& geo) const
    {
      double volume = geo.volume();
      
      typedef typename GeometryType :: ctype coordType;
      enum { dim = GridType :: dimension };
      const ReferenceElement< coordType, dim > & refElem =
             ReferenceElements< coordType, dim >::general(geo.type());
             
      double volRef = refElem.volume();

      double massVolinv = volRef/volume;
      return massVolinv;
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
    const LocalIdSetType & localIdSet_;
    const DiscreteGradientSpaceType & gradientSpace_;
    mutable SingleCommunicationManagerType singleCommunicate_;
    mutable GradCommunicationManagerType gradCommunicate_;
    
    mutable GradDestinationType gradTmp_; 
    mutable GradDestinationType gradRhs_; 
    mutable GradDestinationType * massTmp_; 
    //mutable DestinationType multTmp_;
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


    const int elemOrder_,faceOrder_;
    int volumeQuadOrd_,faceQuadOrd_;

    mutable MatrixHandlerType matrixHandler_;
    // marker for new entities 
    mutable  EntityMarkerType entityMarker_;

    mutable bool matrixAssembled_;
    const bool verbose_;
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

    //! type of restrict and prolong operator during adaptation 
    typedef FEOperatorType RestrictProlongOperatorType;

  private:  
    DiscreteModelType& problem_; 
    DiscreteFunctionSpaceType& spc_;
    
    //ElliptPrevPassType feStartPass_;
    mutable FEOperatorType op_;

    // define type of inverse operator 
    typedef typename DiscreteModelType :: Traits :: Traits :: template 
      InverseOperator<DestinationType,FEOperatorType>:: InverseOperatorType InverseOperatorType;

    typedef typename PreviousPassType :: DestinationType GradDestinationType;

    const double eps_;
    const int maxIterFactor_; 
    mutable int maxIter_;

    InverseOperatorType invOp_; 

    mutable DestinationType rhs_;

    mutable int sequence_;
    const bool verbose_;
      
    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
    mutable CommunicationManagerType comm_;
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param eps epsilon for interative solver 
    //! \param maxIterFactor factor for number of max iterations 
    //! \param verbose if true some output is given 
    LocalDGElliptPass(DiscreteModelType& problem, 
                PreviousPassImp & pass, 
                DiscreteFunctionSpaceType& spc,
                double eps = 1e-10 , int maxIterFactor = 3 , bool verbose = false )
      : BaseType(pass.previousPass(),spc)
      , problem_(problem)
      , spc_(spc) 
      , op_(problem,pass,pass.previousPass(),spc,verbose)
      , eps_(eps)
      , maxIterFactor_(maxIterFactor) 
      , maxIter_( maxIterFactor_ * spc_.size() )
      , invOp_(op_,eps,eps,maxIter_,verbose)
      , rhs_("FEPass::RHS",spc)
      , sequence_(-1)
      , verbose_(verbose)
      , comm_(spc_)
    {
    }

    //! do nothing here 
    void applyLocal(EntityType& en) const
    {
    }
    
    GradDestinationType & tmpMemory() { return op_.tmpMemory(); }

    //! return restrict and prolong operator for fe-pass 
    RestrictProlongOperatorType & restrictProlongOperator () { return op_; }

    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      // for first time buildMatrix 
      if( sequence_ < 0 ) 
      {
        dest.clear();
        op_.prepare(arg,rhs_);
        op_.buildMatrix( arg, rhs_ );
        sequence_ = spc_.sequence();
      }
      else if(sequence_ != spc_.sequence())
      {
        // only clear destination if matrix has really changed 
        // otherwise keep old value as initial value 
        dest.clear();
        op_.prepare(arg,rhs_);

        // for unstructured grids we can use the re-build method 
        if(Capabilities::IsUnstructured<GridType>::v)
        {
          //op_.reBuildMatrix( arg, rhs_ , true );
          op_.buildMatrix( arg, rhs_ );
        }
        else 
        {
          // otherwise just make all new 
          op_.buildMatrix( arg, rhs_ );
        }
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

      // solve the system 
      invOp_(rhs_,dest);

      // do data exchange 
      comm_.exchange ( dest );
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

    typedef CommunicationManager<DiscreteFunctionSpaceType> CommunicationManagerType; 
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
      prevPass_(pass),
      comm_(spc_)
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

    //! calls evalGradient of previous pass 
    void operator () (const GlobalArgumentType& arg, DestinationType& dest) const 
    {
      // normal call procedure 
      prevPass_.pass(arg);

      // now get gradient from previous pass 
      prevPass_.evalGradient(prevPass_.destination(),dest);

      // exchange data 
      comm_.exchange( dest );
    }

    //! nothing to do here
    void applyLocal(EntityType& en) const
    {
    }
    
    //! nothing to do here
    void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
    }

    //! nothing to do here 
    void finalize(const ArgumentType& arg, DestinationType& dest) const
    {
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    DiscreteFunctionSpaceType& spc_;
    mutable PreviousPassImp & prevPass_;
    mutable CommunicationManagerType comm_;
  };

} // end namespace Dune

#endif
