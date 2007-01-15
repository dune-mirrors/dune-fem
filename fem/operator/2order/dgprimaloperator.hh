#ifndef DUNE_DGPRIMALOPERATOR_HH
#define DUNE_DGPRIMALOPERATOR_HH

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
    
    //! type of underlying matrix implementation 
    typedef MatrixObjectImp MatrixObjectType; 

    typedef typename MatrixObjectType::LocalMatrixType LocalMatrixType;
    typedef typename MatrixObjectType::MatrixType MatrixType;
    typedef typename MatrixObjectType::PreconditionMatrixType PreconditionMatrixType;
    
    typedef typename DiscreteModelType :: BoundaryIdentifierType BoundaryIdentifierType;    

    typedef typename LocalIdSetType :: IdType LocalIdType;
    typedef std::set< LocalIdType > EntityMarkerType; 

    typedef GradJacobianRangeType FluxRangeType; 
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
      gradPass_(gradPass),
      gradProblem_(gradPass_.problem()),
      arg_(0),
      dest_(0),
      spc_(spc),
      gridPart_(spc_.gridPart()),
      localIdSet_(gridPart_.grid().localIdSet()),
      gradientSpace_(gradPass_.space()),
      dtMin_(std::numeric_limits<double>::max()),
      time_(0),
      elemOrder_(std::max(spc_.order(),gradientSpace_.order())),
      faceOrder_(std::max(spc_.order(),gradientSpace_.order())+1),
      volumeQuadOrd_(2*elemOrder_),
      faceQuadOrd_(2*faceOrder_),
      matrixObj_(spc_,spc_, problem_.preconditioning()),
      entityMarker_(),
      coeffEn_(1.0),
      coeffNb_(1.0),
      matrixAssembled_(false),
      beta_(0.0),
      bilinearPlus_(true),
      power_( 2*spc_.order() ),
      notBabuskaZlamal_(true),
      betaNotZero_(false)
    {
      if( ! (spc_.order() > 0))
      {
        std::cerr << "DG Primal operator only working for spaces with polynomial order > 0! \n";
        assert(false);
        abort();
      }
      
      // if parameter file is not empty read parameter 
      if(paramFile != "")
      {
        readParameter(paramFile,"beta",beta_);
        int bplus = 1;
        readParameter(paramFile,"B_{+,-}",bplus);
        assert( (bplus == 0) || (bplus == 1) ); 
        bilinearPlus_ = (bplus == 0) ? false : true; 

        int zlamal = 0;
        readParameter(paramFile,"Babuska-Zlamal",zlamal);
        notBabuskaZlamal_ = (zlamal == 1) ? false : true;
      }

      betaNotZero_ = (std::abs(beta_) > 0.0);

      if( !betaNotZero_ && !notBabuskaZlamal_)
      {
        std::cerr << "ERROR: beta == 0.0 and Babuska-Zlamal == 1 !";
        std::cerr << " Choose either beta > 0 or Babuska-Zlamal = 1" << std::endl;
        assert(false);
        exit(1);
      }

      if( !betaNotZero_ )
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
            std::cout << "DGPrimalOperator: using NIPG method, beta = " << beta_ << " !\n"; 
          }
          else 
          {
            std::cout << "DGPrimalOperator: using Interior Penalty method, beta = " << beta_ << " !\n"; 
          }
        }
        else 
        {
          std::cout << "DGPrimalOperator: using Babuska-Zlamal method, beta = " << beta_ << " !\n"; 
        }
      }
      
      //assert( matrixObj_.hasMassMatrix() == gradProblem_.hasSource() );
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

    //! Estimate for the timestep size
    double timeStepEstimate() const {
      return dtMin_;
    }

    //! setup matrix 
    void buildMatrix(const ArgumentType & arg, DestinationType & rhs)
    {
      // reserve memory and clear matrices 
      matrixObj_.reserve(false);

      rhs.clear();
      dest_ = &rhs;

      // build matrix and rhs 
      this->compute( arg, rhs );
      matrixAssembled_ = true;

      // create pre-condition matrix if activated 
      matrixObj_.createPreconditionMatrix();
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
        LocalMatrixType matrixEn(matrixObj_.matrix(),
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
      
      LocalMatrixType matrixEn(matrixObj_.matrix(),
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

      // resize caches 
      if( (int) tau_.size() != numDofs) tau_.resize(numDofs);
      if( (int) phi_.size() != numDofs) phi_.resize(numDofs);
      
      /////////////////////////////////
      // Volumetric integral part
      /////////////////////////////////
      const int quadNop = volQuad.nop();

      RangeType rhsval(0.0);

      // set default value to fMat 
      coeffEn_ = 1.0;
      coeffNb_ = 1.0;
      
      for (int l = 0; l < quadNop ; ++l) 
      {
        // calc integration element 
        const double intel = volQuad.weight(l)
            *geo.integrationElement(volQuad.point(l)) * massVolElInv;

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
        if(problem_.hasCoefficient())
        {
          // call anayltical flux of discrete model 
          caller_.evaluateCoefficient(en, volQuad, l, coeffEn_ );
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

          // apply coefficient 
          if(problem_.hasCoefficient())
          {
            for(int i=0; i<dimRange; ++i)
            {
              psitmp[i] = psi[i];
              psi[i] = 0.0;
              coeffEn_.umv(psitmp[i],psi[i]);
            }
          }

          // add diagonal entry
          {
            double val = bsetEn.evaluateGradientSingle(k, en, volQuad, l, psi )
                         * intel;
            matrixEn.add( k , k , val );
          }

          // add secondary diagonal
          // assume matrix is symectric
          // entry (k,j) == entry (j,k)
          for (int j = k+1; j < numDofs; ++j) 
          {
            // eval grad tau  
            double val = bsetEn.evaluateGradientSingle(j, en, volQuad, l, psi )
                       * intel;
            // add k,j 
            matrixEn.add( k , j , val );

            // add j,k
            matrixEn.add( j , k , val );
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
              problem_.boundaryValue(nit,t,
                faceQuadInner.localPoint(l),boundaryValue);

            ///////////////////////////////
            //  evaluate coefficients 
            ///////////////////////////////
            JacobianRangeType norm(0.0);
            if(problem_.hasCoefficient())
            {
              // call anayltical flux of discrete model 
              caller_.evaluateCoefficientBoundary(nit,faceQuadInner,l,coeffEn_);
              for(int i=0; i<dimRange; ++i)
              {
                coeffEn_.umv(unitNormal,norm[i]);
              }
            }
            else 
            {
              for(int i=0; i<dimRange; ++i)
              {
                norm[i] = unitNormal;
              }
            }

            // cache base functions evaluations
            for(int k=0; k<numDofs; ++k)
            { 
              tau_[k] = bsetEn.evaluateGradientSingle(k,en,faceQuadInner,l, norm);  
              bsetEn.eval(k,faceQuadInner,l, phi_[k]);
            }
               
            for(int k=0; k<numDofs; ++k)
            { 
              // if not Babuska-Zlamal method, add boundary terms 
              if( notBabuskaZlamal_ )
              {
                if( bndType.isDirichletNonZero())
                {
                  // only valid for dim range = 1
                  double rhsVal1 = boundaryValue[0] * tau_[k];

                  rhsVal1 *= bilinIntel;
                  singleRhs[k] += rhsVal1;
                }

                // grad w * v 
                // only on non Neumann type boundaries
                if( ! bndType.isNeumannType() )
                {
                  for (int j = 0; j < numDofs; ++j) 
                  {
                    double val = tau_[j] * phi_[k][0];
                    val *= -intel;
                    matrixEn.add( k , j , val );
                  }
                }
                  
                // w * grad v
                if(bndType.isDirichletType())
                {
                  for (int j = 0; j < numDofs; ++j) 
                  {
                    double val = tau_[k] * phi_[j][0];
                    val *= bilinIntel;
                    matrixEn.add( k , j , val );
                  }
                }
              }
                
              // dirichlet boundary values for u 
              if(bndType.isNeumannNonZero())
              {
                // only valid for dim range = 1
                double rhsVal1 = boundaryValue[0] * phi_[k][0];

                rhsVal1 *= intelFactor;
                singleRhs[k] += rhsVal1;
              }

              if( betaNotZero_ )
              {
                // stabilization 
                if( bndType.isDirichletType())
                {
                  for (int j = 0; j < numDofs; ++j) 
                  {
                    double phiVal = phi_[k] * phi_[j]; 
                    phiVal *= facBeta;
                    matrixEn.add( k , j , phiVal );
                  }
                }
                
                // dirichlet boundary values for u 
                if(bndType.isDirichletNonZero())
                {
                  // only valid for dim range = 1
                  double rhsVal1 = boundaryValue[0] * phi_[k];
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
                            EntityType & en, const double massVolElInv,
                            VolumeQuadratureType & volQuad,
                            const QuadratureImp & faceQuadInner, 
                            const QuadratureImp & faceQuadOuter, 
                            const BaseFunctionSetType & bsetEn, 
                            LocalMatrixType & matrixEn) const
    {
      const int numDofs = bsetEn.numBaseFunctions();

      EntityPointerType neighEp = nit.outside();
      EntityType&            nb = *neighEp;

      caller_.setNeighbor(nb);

      // create matrix handles for neighbor 
      LocalMatrixType matrixNb(matrixObj_.matrix(),
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

        RangeType phi_j;
        RangeType phiNeigh;

        JacobianRangeType normEn(0.0);
        JacobianRangeType normNb(0.0);
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
            coeffEn_.umv(unitNormal,normEn[i]);
            coeffNb_.umv(unitNormal,normNb[i]);
          }
        }
        else 
        {
          for(int i=0; i<dimRange; ++i)
          {
            normEn[i] = unitNormal;
            normNb[i] = unitNormal;
          }
        }
               
        // cache base functions evaluations
        for(int k=0; k<numDofs; ++k)
        { 
          // eval base functions 
          bsetEn.eval(k,faceQuadInner,l, phi_[k]);
          // eval gradient 
          tau_[k] = bsetEn.evaluateGradientSingle(k,en,faceQuadInner,l, normEn);  
        }
               
        for(int k=0; k<numDofs; ++k)
        {
          // this terms dissapear if Babuska-Zlamal is used 
          if(notBabuskaZlamal_)
          {
            for(int j=0; j<numDofs; ++j)
            {
              double tauNeigh = bsetNeigh.evaluateGradientSingle(j,nb,faceQuadOuter,l, normNb);      

              // v^+ * (grad w^+  + grad w^-)
              {
                numericalFlux2(phi_[k] , tau_[j] , tauNeigh, resultLeft, resultRight);

                double valLeft = resultLeft[0];
                valLeft *= -intel;

                matrixEn.add( k , j , valLeft );

                double valRight = resultRight[0];
                valRight *= -intel;

                matrixNb.add( k , j , valRight );
              }

              bsetNeigh.eval(j,faceQuadOuter,l, phiNeigh );      
              
              // v^+ * (grad w^+  + grad w^-)
              {
                numericalFlux(tau_[k] , phi_[j] , phiNeigh , resultLeft, resultRight);

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
              phi_j = bsetEn.evaluateSingle(j,faceQuadInner,l, phi_[k] );      
              phiNeigh = bsetNeigh.evaluateSingle(j,faceQuadOuter,l, phi_[k] );      
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

      MatrixAddHandleType massMatrixEn (matrixObj_.massMatrix(),
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
    DGPrimalOperator();
    //! copy constructor not defined 
    DGPrimalOperator(const DGPrimalOperator&);

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
    
    mutable double dtMin_;
  
    TimeProvider* time_;

    const int elemOrder_,faceOrder_;
    int volumeQuadOrd_,faceQuadOrd_;

    mutable MatrixObjectType matrixObj_;
    // marker for new entities 
    mutable  EntityMarkerType entityMarker_;

    // return type of analyticalFlux 
    mutable FluxRangeType coeffEn_;
    mutable FluxRangeType coeffNb_;
    // caches for base function evaluation 
    mutable std::vector<RangeFieldType> tau_;
    mutable std::vector<RangeType> phi_;

    mutable bool matrixAssembled_;
    double beta_;

    // if true B_+ is used otherwise B_-
    bool bilinearPlus_;
    double power_;
    bool notBabuskaZlamal_;
    bool betaNotZero_;
  };
  
} // end namespace Dune
#endif
