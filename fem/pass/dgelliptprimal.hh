#ifndef DUNE_NEWELLIPTPASS_HH
#define DUNE_NEWELLIPTPASS_HH

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
  ////////////////////////////////////////////////////////////
  //
  //  --LocalDGPrimalOperator 
  //
  ////////////////////////////////////////////////////////////
  //! Concrete implementation of Pass for LDG.
  template <class DiscreteModelImp, class GradientPassImp, 
            class PreviousPassImp>
  class LocalDGPrimalOperator 
    : public LocalPass<DiscreteModelImp, PreviousPassImp> 
    , public OEMSolver::PreconditionInterface
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass<DiscreteModelImp, PreviousPassImp> BaseType;

    typedef LocalDGPrimalOperator<DiscreteModelImp,GradientPassImp,PreviousPassImp> ThisType;

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
    typedef DiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::DimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::DimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    enum { dim = GridType :: dimension };
    

    typedef FieldMatrix<double,rows,rows> TensorType;
    typedef FieldMatrix<double,dim,dim> JacobianInverseType;
    
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
  public:
    //- Public methods
    //! Constructor
    //! \param problem Actual problem definition (see problem.hh)
    //! \param pass Previous pass
    //! \param spc Space belonging to the discrete function local to this pass
    //! \param volumeQuadOrd defines the order of the volume quadrature which is by default 2* space polynomial order 
    //! \param faceQuadOrd defines the order of the face quadrature which is by default 2* space polynomial order 
    LocalDGPrimalOperator(DiscreteModelType& problem, 
                GradientPassType & gradPass,
                PreviousPassType& pass, 
                const DiscreteFunctionSpaceType& spc,
                bool verbose = false, 
                int volumeQuadOrd =-1,int faceQuadOrd=-1) 
      : BaseType(pass, spc),
      caller_(problem),
      problem_(problem),
      gradPass_(gradPass),
      gradProblem_(gradPass_.problem()),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      localIdSet_(gridPart_.grid().localIdSet()),
      gradientSpace_(gradPass_.space()),
      singleCommunicate_(spc_),
      dtMin_(std::numeric_limits<double>::max()),
      time_(0),
      elemOrder_(std::max(spc_.order(),gradientSpace_.order())),
      faceOrder_(std::max(spc_.order(),gradientSpace_.order())+1),
      volumeQuadOrd_( (volumeQuadOrd < 0) ? (2*elemOrder_) : volumeQuadOrd ),
      faceQuadOrd_( (faceQuadOrd < 0) ? (2*faceOrder_) : faceQuadOrd ),
      matrixHandler_(spc_,gradientSpace_,
                     gradProblem_.hasSource(),problem_.preconditioning()),
      entityMarker_(),
      matrixAssembled_(false),
      verbose_(verbose),
      beta_(0.0),
      bilinearPlus_(true),
      power_( 2*spc_.order() ),
      zlamal_(false),
      betaNotZero_(false)
    {
      readParameter("parameter","beta",beta_);
      int bplus = 1;
      readParameter("parameter","B{+,-}",bplus);
      assert( (bplus == 0) || (bplus == 1) ); 
      bilinearPlus_ = (bplus == 0) ? false : true; 

      int zlamal = 0;
      readParameter("parameter","zlamal",zlamal);
      zlamal_ = (zlamal == 1) ? true : false;

      betaNotZero_ = (std::abs(beta_) > 0.0);
      assert( (zlamal_) ? (std::abs(beta_) > 0.0) : 1);
      
      assert( matrixHandler_.hasMassMatrix() == gradProblem_.hasSource() );
      assert( volumeQuadOrd_ >= 0 );
      assert( faceQuadOrd_ >= 0 );

      if( ! (spc_.order() > 0))
      {
        std::cerr << "DG Primal operator only working for spaces with polynomial order > 0! \n";
        assert(false);
        abort();
      }
      
      if(problem_.hasSource())
      {
        std::cerr << "Source for DGElliptPass not supported yet. \n";
        abort();
      }
    }

    //! Destructor
    virtual ~LocalDGPrimalOperator() 
    {
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
      // reserve memory and clear matrices 
      matrixHandler_.reserve(verbose_);

      rhs.clear();
      dest_ = &rhs;

      // build matrix and rhs 
      this->compute( arg, rhs );
      matrixAssembled_ = true;

      // create pre-condition matrix if activated 
      createPreconditionMatrix();
    }

    //! do matrix vector multiplication, used by InverseOp  
    void operator () (const DestinationType & arg, DestinationType& dest) const 
    {
      multOEM( arg.leakPointer(), dest.leakPointer());
    }
    
    //! do matrix vector multiplication, used by OEM-Solver and DuneODE Solvers  
    // --multOEM
    void multOEM(const double * arg, double * dest) const
    {
      // calc dest = stabMatrix * arg 
      matrixHandler_.stabMatrix().multOEM(arg,dest);
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
    const PreconditionMatrixType & preconditionMatrix () const { 
      return matrixHandler_.pcMatrix(); 
    }

    //! returns true if preconditioning matrix has been build 
    bool hasPreconditionMatrix() const  { 
      return matrixHandler_.hasPcMatrix(); 
    }
           
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      caller_.setArgument(*arg_);
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
        MatrixAddHandleType matrixEn(matrixHandler_.stabMatrix(),
                                       en, spc_, en, spc_ ); 
        matrixEn.clear();
      }
    }
    
    //! apply operator on entity 
    void applyLocal(EntityType& en) const
    {
      // local function for right hand side 
      typedef typename DestinationType :: LocalFunctionType SingleLFType; 
      SingleLFType singleRhs = dest_->localFunction(en); //rhs
      
      MatrixAddHandleType matrixEn(matrixHandler_.stabMatrix(),
                                     en, spc_, en, spc_ ); 
      
      // make entities known in callers
      caller_.setEntity(en);

      VolumeQuadratureType volQuad(en, volumeQuadOrd_);

      const GeometryType & geo = en.geometry();
      //const double volume = geo.volume();
      const double massVolElInv = massVolumeInv(geo);

      // get base function set of single space 
      const BaseFunctionSetType& bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////
      const int quadNop = volQuad.nop();

      RangeType rhsval(0.0);
      
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc factor for bas functions 
        const double rhsFactor = volQuad.weight(l)
            *geo.integrationElement(volQuad.point(l));

        const double intel = rhsFactor * massVolElInv;

        const JacobianInverseType& inv =
          geo.jacobianInverseTransposed(volQuad.point(l));

        ////////////////////////////////////
        // create rightHandSide
        ////////////////////////////////////
        {
          // set default value 

          // eval rightHandSide function 
          // if empty, rhs stays 0.0
          caller_.rightHandSide(en, volQuad, l, rhsval );

          // scale with quadrature weight  
          rhsval *= volQuad.weight(l); 

          for (int j = 0; j < numDofs; ++j) 
          {
            singleRhs[j] += bsetEn.evaluateSingle(j, volQuad, l, rhsval );
          }
        }
        
        ///////////////////////////////
        //  evaluate coefficients 
        ///////////////////////////////
        if(problem_.hasFlux())
        {
          // call anayltical flux of discrete model 
          //caller_.analyticalFlux(en, volQuad, l, fMat_ );
        }

        JacobianRangeType psitmp(0.0);
        for(int k = 0; k < numDofs; ++k)
        {
          // eval grad psi 
          bsetEn.jacobian( k, volQuad, l, psitmp );

          JacobianRangeType psi(0.0);
          // apply inverse jacobian 
          for(int i=0; i<dimRange; ++i) 
          {
            inv.umv(psitmp[i], psi[i]);
          }

          for (int j = 0; j < numDofs; ++j) 
          {
            // eval grad tau  
            double val = bsetEn.evaluateGradientSingle(j, en, volQuad, l, psi )
                       * intel;
            matrixEn.add( k , j , val );
          }
        }
      } // end element integral 
  
      /////////////////////////////////
      // Surface integral part
      /////////////////////////////////
      
      IntersectionIteratorType endnit = gridPart_.iend(en); 
      for (IntersectionIteratorType nit = gridPart_.ibegin(en); nit != endnit; ++nit) 
      { 
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
                  bsetEn,matrixEn);
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
                  bsetEn,matrixEn);
          }
        } // end if neighbor 

        // if intersection with boundary 
        if (nit.boundary()) 
        { 
          FaceQuadratureType faceQuadInner(gridPart_, nit, faceQuadOrd_,
                                           FaceQuadratureType::INSIDE);

          // get time if time ptovider exists  
          const double t = (time_) ? (time_->time()) : 0.0;

          RangeType phi_j(0.0);
          RangeType phi_k(0.0);

          const int quadNop = faceQuadInner.nop();
          for (int l = 0; l < quadNop ; ++l) 
          {
            DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
            const double faceVol = unitNormal.two_norm();
            unitNormal *= 1.0/faceVol;

            const double intelFactor = faceQuadInner.weight(l) * massVolElInv;
            
            // integration element factor 
            const double intel = intelFactor * faceVol;

            // intel switching between bilinear from B_+ and B_-  
            const double bilinIntel = (bilinearPlus_) ? intel : -intel;

            // overall beta factor 
            const double facBeta = factorBeta(intelFactor,faceVol);

            // get boundary value 
            RangeType boundaryValue(0.0);

            BoundaryIdentifierType bndType = 
              problem_.boundaryValue(nit,unitNormal,t,
                faceQuadInner.localPoint(l),boundaryValue);

            JacobianRangeType norm;
            for(int i=0; i<dimRange; ++i)
            {
              norm[i] = unitNormal;
            }
               
            {
              for(int k=0; k<numDofs; ++k)
              { 
                bsetEn.eval(k,faceQuadInner,l, phi_k);
                double tau_k = bsetEn.evaluateGradientSingle(k,en,faceQuadInner,l, norm);

                // if not Babuska-Zlamal method, add boundary terms 
                if( !zlamal_ )
                {
                  if( bndType.isDirichletNonZero())
                  {
                    // only valid for dim range = 1
                    double rhsVal1 = boundaryValue[0] * tau_k;

                    rhsVal1 *= bilinIntel;
                    singleRhs[k] += rhsVal1;
                  }

                  // dirichlet boundary values for u 
                  if(bndType.isNeumannNonZero())
                  {
                    //std::cout << "Adding neumann value " << boundaryValue[0] << "\n";
                    // only valid for dim range = 1
                    double rhsVal1 = boundaryValue[0] * phi_k[0];

                    rhsVal1 *= intel;
                    singleRhs[k] -= rhsVal1;
                  }

                  for (int j = 0; j < numDofs; ++j) 
                  {
                    bsetEn.eval(j, faceQuadInner, l, phi_j );
                    double tau_j = bsetEn.evaluateGradientSingle(j,en,faceQuadInner,l, norm);
                  
                    // grad v * w 
                    {
                      double val = tau_j * phi_k[0];
                      val *= -intel;
                      matrixEn.add( k , j , val );
                    }
                    
                    // w * grad v
                    {
                      double val = tau_k * phi_j[0];
                      val *= bilinIntel;
                      matrixEn.add( k , j , val );
                    }
                  }
                }
                  
                if( betaNotZero_ )
                {
                  // stabilization 
                  for (int j = 0; j < numDofs; ++j) 
                  {
                    double phiVal = bsetEn.evaluateSingle(j, faceQuadInner, l, phi_k);
                    phiVal *= facBeta;
                    matrixEn.add( k , j , phiVal );
                  }
                  
                  // dirichlet boundary values for u 
                  if(bndType.isDirichletNonZero())
                  {
                    // only valid for dim range = 1
                    double rhsVal1 = boundaryValue[0] * phi_k;
                    rhsVal1 *= facBeta;
                    singleRhs[k] += rhsVal1;
                  }
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
      if(zlamal_)
      {
        return (beta_ * intelFactor * faceVol * pow(faceVol , -power_ )); 
      }
      else 
      {
        return (beta_ * intelFactor);
      }
    }

    template <class QuadratureImp> 
    void applyLocalNeighbor(IntersectionIteratorType & nit, 
                            EntityType & en, const double massVolElInv,
                            VolumeQuadratureType & volQuad,
                            const QuadratureImp & faceQuadInner, 
                            const QuadratureImp & faceQuadOuter, 
                            const BaseFunctionSetType & bsetEn, 
                            MatrixAddHandleType & matrixEn) const
    {
      const int numDofs = bsetEn.numBaseFunctions();

      EntityPointerType neighEp = nit.outside();
      EntityType&            nb = *neighEp;

      caller_.setNeighbor(nb);

      // create matrix handles for neighbor 
      MatrixAddHandleType matrixNb(matrixHandler_.stabMatrix(),
                                     en, spc_, nb, spc_ ); 
      
      // get base function set 
      const BaseFunctionSetType& bsetNeigh = spc_.baseFunctionSet(nb);
     
      const int quadNop = faceQuadInner.nop();
      for (int l = 0; l < quadNop ; ++l) 
      {
        DomainType unitNormal(nit.integrationOuterNormal(faceQuadInner.localPoint(l)));
        const double faceVol = unitNormal.two_norm();
        unitNormal *= 1.0/faceVol; 

        // integration element factor 
        const double intelFactor = faceQuadInner.weight(l) * massVolElInv; 
        const double intel = faceVol * intelFactor; 

        // intel switching between bilinear from B_+ and B_-  
        const double bilinIntel = (bilinearPlus_) ? intel : -intel;

        // overall beta factor 
        const double facBeta = factorBeta(intelFactor,faceVol);

        ////////////////////////////////////////////////////////////
        RangeType resultLeft(0.0);
        RangeType resultRight(0.0);

        RangeType phi_k;
        RangeType phi_j;
        RangeType phiNeigh;

        JacobianRangeType tau2;
        JacobianRangeType norm; 
        for(int i=0; i<dimRange; ++i)
        {
          norm[i] = unitNormal;
        }
               
        for(int k=0; k<numDofs; ++k)
        {
          bsetEn.eval(k,faceQuadInner,l, phi_k );      
          double tau_k = bsetEn.evaluateGradientSingle(k, en, faceQuadInner,l, norm );
          
          // this terms dissapear if Babuska-Zlamal is used 
          if(!zlamal_)
          {
            for(int j=0; j<numDofs; ++j)
            {
              double tau_j = bsetEn.evaluateGradientSingle(j,en,faceQuadInner,l, norm);      
              double tauNeigh = bsetNeigh.evaluateGradientSingle(j,nb,faceQuadOuter,l, norm);      

              // v^+ * (grad w^+  + grad w^-)
              {
                numericalFlux2(phi_k , tau_j , tauNeigh, resultLeft, resultRight);

                double valLeft = resultLeft[0];
                valLeft *= -intel;

                matrixEn.add( k , j , valLeft );

                double valRight = resultRight[0];
                valRight *= -intel;

                matrixNb.add( k , j , valRight );
              }

              bsetEn.eval(j,faceQuadInner,l, phi_j );      
              bsetNeigh.eval(j,faceQuadOuter,l, phiNeigh );      
              
              // v^+ * (grad w^+  + grad w^-)
              {
                numericalFlux(tau_k , phi_j , phiNeigh , resultLeft, resultRight);

                double valLeft = resultLeft;
                valLeft *= bilinIntel;

                matrixEn.add( k , j , valLeft );

                double valRight = resultRight;
                valRight *= bilinIntel;

                matrixNb.add( k , j , valRight );
              }
            }
          }

          if( betaNotZero_ )
          {
            for(int j=0; j<numDofs; ++j)
            {
              phi_j = bsetEn.evaluateSingle(j,faceQuadInner,l, phi_k );      
              phiNeigh = bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phi_k );      
              {
                RangeType resLeft, resRight;
                
                numericalFluxStab(phi_j, phiNeigh , resLeft, resRight);

                double valLeft = facBeta;
                valLeft  *= resLeft;

                matrixEn.add( k , j , valLeft );

                double valRight = facBeta;
                valRight *= resRight;

                matrixNb.add( k , j , valRight );
              }
            }
          }
        }
      }
    }

    void numericalFlux(const double & grad, 
                       const RangeType & phiLeft,
                       const RangeType & phiRight, 
                       RangeType & resultLeft,
                       RangeType & resultRight) const
    {
      resultLeft  = grad;
      resultLeft *= 0.5;
      
      resultRight  = grad;
      resultRight *= 0.5; 

      resultLeft  *= phiLeft;
      resultRight *= -phiRight;
    }
                       
    void numericalFlux2(const RangeType & phi,
                        const double & gradLeft, 
                        const double & gradRight,
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
      resultLeft  =  phiLeft; 
      resultRight = -phiRight; 
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
        PreconditionMatrixType & diag = matrixHandler_.pcMatrix(); 
        diag.clear();
        
        /*
        if(gradProblem_.hasSource())
        {
          matrixHandler_.divMatrix().getDiag( matrixHandler_.massMatrix(), matrixHandler_.gradMatrix() , diag );
        }
        else 
        {
          matrixHandler_.divMatrix().getDiag( matrixHandler_.gradMatrix() , diag );
        }
        */

        matrixHandler_.stabMatrix().addDiag( diag );

        double * diagPtr = diag.leakPointer();
        const int singleSize = spc_.size();
        for(register int i=0; i<singleSize; ++i) 
        {
          double val = diagPtr[i]; 
          // when using parallel Version , we could have zero on diagonal
          // for ghost elements 
          assert( (spc_.grid().comm().size() > 1) ? 1 : (std::abs( val ) > 0.0 ) );
          if( std::abs( val ) > 0.0 )
          {
            val = 1.0/val; 
            diagPtr[i] = val;
          }
        }
      }
    }

    void updateLocal(EntityType& en) const
    {
      /*
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
      */
    }

  private:
    // needs to be friend for conversion check 
    friend class Conversion<ThisType,OEMSolver::PreconditionInterface>;
    //! empty constructor not defined 
    LocalDGPrimalOperator();
    //! copy constructor not defined 
    LocalDGPrimalOperator(const LocalDGPrimalOperator&);

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
      //return 1.0;
    }

  private:
    mutable DiscreteModelCallerType caller_;
    DiscreteModelType& problem_; 
    GradientPassType & gradPass_; 
    GradientDiscreteModelType& gradProblem_; 
    
    mutable ArgumentType* arg_;
    mutable DestinationType* dest_;

    const DiscreteFunctionSpaceType& spc_;
    const GridPartType & gridPart_;
    const LocalIdSetType & localIdSet_;
    const DiscreteGradientSpaceType & gradientSpace_;
    mutable SingleCommunicationManagerType singleCommunicate_;
    
    mutable double dtMin_;
  
    TimeProvider* time_;

    const int elemOrder_,faceOrder_;
    int volumeQuadOrd_,faceQuadOrd_;

    mutable MatrixHandlerType matrixHandler_;
    // marker for new entities 
    mutable  EntityMarkerType entityMarker_;

    mutable bool matrixAssembled_;
    const bool verbose_;
    double beta_;

    // if true B_+ is used otherwise B_-
    bool bilinearPlus_;
    double power_;
    bool zlamal_;
    bool betaNotZero_;
  };
  
} // end namespace Dune
#endif
