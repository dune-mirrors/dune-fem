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

#include <dune/fem/misc/boundaryidentifier.hh>
#include <dune/fem/solver/timeprovider.hh>
#include <dune/fem/solver/oemsolver/preconditioning.hh>

#include <dune/fem/space/combinedspace.hh>

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
**   div(A(x)\nabla u) &=& f(x)  \quad\mbox{in}\quad \Omega    \\
** \f}
** where \f$ v \f$ is to be computed.
** @{
**************************************************************************/
  ////////////////////////////////////////////////////////////
  //
  //  --DGPrimalOperator 
  //
  ////////////////////////////////////////////////////////////
  /** \brief Operator assembling matrix for DG methods for elliptic problems. 
      Currently implemented are:
        - Interior Penalty
        - Baumann-Oden
        - NIPG 
        - Babuska-Zlamal
        - Compact LDG (CDG) 

        References:
          The first 4 methods can be found for example in:
            D.N. Arnold, F. Brezzi, B. Cockburn, L.D. Marini: Unified
            Analysis of Discontinuous Galerkin Methods for Elliptic
            Problems SIAM J. Num. Anal, 39  (2002), 1749-1779.
            http://www.imati.cnr.it/~marini/reports/dgm_siam.ps.gz

          The Compact LDG method is described in detail in:
            J. Peraire and P.-O. Persson,
            The Compact Discontinuous Galerkin (CDG) Method for Elliptic Problems.
            SIAM J. Sci. Comput., to appear.
            http://www.mit.edu/~persson/pub/peraire07cdg.pdf
  */
  template <class DiscreteModelImp, class GradientPassImp, 
            class PreviousPassImp, class MatrixObjectTraits,
            int pId = -1 >
  class DGPrimalOperator 
    : public LocalPass< DiscreteModelImp, PreviousPassImp, pId > 
    , public OEMSolver::PreconditionInterface
  {
  public:
    //- Typedefs and enums
    //! Base class
    typedef LocalPass< DiscreteModelImp, PreviousPassImp, pId > BaseType;

    typedef DGPrimalOperator<DiscreteModelImp,GradientPassImp,
            PreviousPassImp,MatrixObjectTraits> ThisType;

    typedef GradientPassImp GradientPassType; 
    typedef typename GradientPassType :: DiscreteModelType
      GradientDiscreteModelType;

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
    typedef typename DiscreteFunctionSpaceType::RangeType RangeType;
    typedef typename DiscreteFunctionSpaceType::JacobianRangeType JacobianRangeType;
    typedef typename DiscreteFunctionSpaceType :: GridType GridType;
        
    // Range of the destination
    enum { dimDomain = DiscreteFunctionSpaceType::dimDomain };
    enum { dimRange = DiscreteFunctionSpaceType::dimRange };
    enum { cols = JacobianRangeType :: cols };
    enum { rows = JacobianRangeType :: rows };
    enum { dim = GridType :: dimension };
    
    enum { dimGradRange = dimDomain * dimRange };
    enum { polOrd = DiscreteFunctionSpaceType::polynomialOrder };

    //! space of gradients of function 
    typedef CombinedSpace< DiscreteFunctionSpaceType, 
                           dimGradRange,
                           PointBased > DiscreteGradientSpaceType; 

    // Types extracted from the discrete function space type
    typedef typename DiscreteFunctionSpaceType::GridPartType GridPartType;
    typedef typename GridPartType :: GridType :: Traits :: LocalIdSet LocalIdSetType; 
    typedef typename DiscreteFunctionSpaceType::DomainType DomainType;
    typedef FieldVector<double,DomainType::dimension-1> FaceDomainType;

    // Types extracted from the underlying grid
    typedef typename GridPartType::IntersectionIteratorType IntersectionIteratorType;
    typedef typename GridType::template Codim<0>::Geometry GeometryType;

    // Various other types
    typedef typename DestinationType::LocalFunctionType LocalFunctionType;
    typedef typename DiscreteModelType::SelectorType SelectorType;
    typedef EllipticDiscreteModelCaller<
      DiscreteModelType, ArgumentType, SelectorType> DiscreteModelCallerType;

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
    typedef TemporaryLocalFunction< DiscreteGradientSpaceType > TemporaryLocalFunctionType;
    typedef MutableArray< std::auto_ptr< TemporaryLocalFunctionType > > TemporaryLocalFunctionArrayType;

    typedef DGMatrixTraits< MatrixObjectTraits > MyOperatorTraits;
    //! type of underlying matrix implementation 
    typedef typename MatrixObjectTraits :: template MatrixObject<
      MyOperatorTraits > :: MatrixObjectType MatrixObjectType;
    //typedef MatrixObjectImp MatrixObjectType; 

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
     \param[in]  gradPass  types and discrete model for the gradient pass
     \param pass Previous pass
     \param spc Space belonging to the discrete function local to this pass
     \param paramFile parameter file to read necessary parameters, if empty 
             default parameters will be applied 
    
     \note Available methods are (chosen by parameters B_{+,-}, beta, and CDG-BZ)
          - Interior Penalty : B_{+,-}: 0 , beta: > 0 (big) , CDG-BZ: 0 
          - Baumann-Oden     : B_{+,-}: 1 , beta: = 0       , CDG-BZ: 0 (needs polOrd > 1) 
          - NIPG             : B_{+,-}: 1 , beta: > 0       , CDG-BZ: 0
          - Babuska-Zlamal   : B_{+,-}: 1 , beta: > 0       , CDG-BZ: 1
          - Compact LDG (CDG): B_{+,-}: 0 , beta: > 0       , CDG-BZ: 1
    */         
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
      gradientSpace_(gridPart_),
      localIdSet_(gridPart_.grid().localIdSet()),
      gridWidth_ ( GridWidthProviderType :: getObject( &spc_.grid())),
      time_(0),
      volumeQuadOrd_(2* spc_.order() ),
      faceQuadOrd_(2*spc_.order() + 2),
      matrixObj_(spc_,spc_, paramFile ),
      coeffEn_(1.0),
      coeffNb_(1.0),
      matrixAssembled_(false),
      betaFactor_(0.0),
      globalBeta_(0.0),
      beta_(0.0),
      bilinearPlus_(true),
      power_( 2*spc_.order() ),
      notBabuskaZlamal_(true),
      compactLDG_(false),
      betaNotZero_(false)
    {
#ifndef DG_DOUBLE_FEATURE
      // we need global upwind vector to select edges 
      upwind_ = M_PI;
      if( dim > 1 ) upwind_[1] = M_LN2;
      if( dim > 2 ) upwind_[2] = M_E;
#endif

      for( int i=0; i<FluxRangeType :: rows ; ++i) 
      {
        // set diagonal to 1 
        coeffEn_[i][i] = coeffNb_[i][i] = 1.0;
        // all others to zero 
        for(int j=i+1; j< FluxRangeType :: cols; ++j) 
        {
          // set default value to fMat 
          coeffEn_[i][j] = coeffNb_[i][j] = 0.0;
          coeffEn_[j][i] = coeffNb_[j][i] = 0.0;
        }
      }

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

        int zlamal = 0;
        if( ! readParameter(paramFile,"CDG-BZ",zlamal, output) )
        {
          if( readParameter(paramFile,"Babuska-Zlamal",zlamal, output) )
          {
            std::cerr << std::endl;
            std::cerr << "WARNING: parameter `Babuska-Zlamal' is deprecated, change it to `CDG-BZ' please in file: ";
            std::cerr << paramFile << " !" << std::endl << std::endl;  
          }
          else 
          {
            success = false; 
          }
        }

        notBabuskaZlamal_ = (zlamal == 1) ? false : true;
      }

      if( ! success )
      {
        if( output ) 
        {
          std::cerr << "\nERROR: Couldn't read parameter! \n";
          std::cerr << "DGPrimalOperator -- Available Methods:\n";
          std::cerr << "Interior Penalty : B_{+,-}: 0 , beta: > 0 (big) , CDG-BZ: 0 \n";
          std::cerr << "Baumann-Oden     : B_{+,-}: 1 , beta: = 0       , CDG-BZ: 0 \n";
          std::cerr << "NIPG             : B_{+,-}: 1 , beta: > 0       , CDG-BZ: 0 \n";
          std::cerr << "Babuska-Zlamal   : B_{+,-}: 1 , beta: > 0       , CDG-BZ: 1 \n";
          std::cerr << "Compact LDG (CDG): B_{+,-}: 0 , beta: > 0       , CDG-BZ: 1 \n\n";
        }
        exit(1);
      }

      betaNotZero_ = (std::abs(betaFactor_) > 0.0);

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
    virtual void setTime ( const double time )
    {
      time_ = time;
    }

    //! Estimate for the timestep size, return value is 0 
    double timeStepEstimate() const {
      return 0.0;
    }

    void compute(const ArgumentType& arg, DestinationType& dest) const
    {
      prepare(arg, dest);

      IteratorType endit = spc_.end();
      for (IteratorType it = spc_.begin(); it != endit; ++it) 
      {
        applyLocal(*it);
      }
            
      // also calculate matrix for overlap 
      if( spc_.grid().overlapSize(0) > 0 ) 
      {
        typedef typename GridPartNewPartitionType<GridPartType,Overlap_Partition>:: NewGridPartType NewGridPartType;
        typedef typename NewGridPartType :: template Codim<0> :: IteratorType IteratorType;

        NewGridPartType gridPart( gridPart_.grid() );
        IteratorType endit = gridPart. template end<0>();
        for(IteratorType it = gridPart. template begin<0>();
            it != endit; ++it)
        {
          applyLocal(*it);
        }
      }

      finalize(arg, dest);
    }

    //! compute matrix entries 
    void computeMatrix(const ArgumentType & arg, DestinationType & rhs)
    {
      // prepare 
      prepare( arg, rhs );

      // if grid has changed, then matrix structure is re-build
      matrixObj_.reserve();

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
      return matrixObj_.preconditionMatrix(); 
    }

    //! returns true if preconditioning matrix has been build 
    bool hasPreconditionMatrix() const  { 
      return matrixObj_.hasPreconditionMatrix(); 
    }
           
    //! In the preparations, store pointers to the actual arguments and 
    //! destinations. Filter out the "right" arguments for this pass.
    virtual void prepare(const ArgumentType& arg, DestinationType& dest) const
    {
      arg_ = const_cast<ArgumentType*>(&arg);
      dest_ = &dest;
      caller_.setArgument(*arg_);

      // only calculate in case of beta not zero 
      if( betaNotZero_ )
      {
        // calculate beta = O(1/h)
        globalBeta_ = betaFactor_/gridWidth_.gridWidth();
      }

      caller_.setTime( time_ );
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

    /////////////////////////////////
    //! apply operator on entity 
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

      // get base function set of single space 
      const BaseFunctionSetType bsetEn = spc_.baseFunctionSet(en);
      const int numDofs = bsetEn.numBaseFunctions();
      assert( numDofs > 0 );

      // set matrix to id matrix 
      if( en.partitionType() != InteriorEntity ) 
      {
        for(int k=0; k<numDofs; ++k) 
        {
          matrixEn.set( k, k, 1);
        }
        return ;
      }

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
          const bool ghostEntity = 
            ( nb.partitionType() != InteriorEntity );
          // only once per intersection or when outside is not interior 
          if( (localIdSet_.id(en) < localIdSet_.id(nb)) 
              || ghostEntity
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
                    //, true
                    , ! ghostEntity 
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

          typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
          const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
          // get number of base functions for gradient space 
          const int numGradBase = enSet.numBaseFunctions();

          // swtich for adding compact LDG lifting operator 
          bool addCompactDG = compactLDG_;
          if( compactLDG_ ) 
          {
            // resize and reset temporary functions 
            resizeTemporaryFunctions( en, r_e_, numDofs, numGradBase );
          } // end compact LDG 

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
              tau_[k] = bsetEn.evaluateGradientSingle(k,en, faceQuadInner[l] , norm);  
              // evaluate phi 
              bsetEn.evaluate(k,faceQuadInner[l] , phi_[k]);
            }

            // only add compact LDG values on dicrichlet boundary 
            addCompactDG = (compactLDG_ && bndType.isDirichletType() );
            if( addCompactDG )
            {
              GradRangeType& tmp = rRets_[0];

              // get numbre of base functions 
              for(int m=0; m<numGradBase; ++m)
              {  
                // eval base functions 
                enSet.evaluate(m, faceQuadInner[l], tmp );
                // apply unit normal 
                eta_[m] = tmp * unitNormal;
              }

              RangeType bndPhi(0.0); 
              // only add this when we have a dirichlet boundary 
              if( bndType.isDirichletType() )
              {
                bndPhi = boundaryValue;
              }

              // calculate coefficients for r_D( g_D - u_h )
              for(int k=0; k<numDofs; ++k)
              {
                RangeType phiVal ( bndPhi );
                phiVal -= phi_[k];
                //RangeType phiVal ( phi_[k] );
                //phiVal -= bndPhi;
                phiVal *= intel;

                TemporaryLocalFunctionType& r_e = *(r_e_[k]);
                // calculate coefficients 
                for(int m=0; m<numGradBase; ++m)
                {
                  r_e[m] -= phiVal * eta_[m]; 
                }
              }
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

          // now add lifting operators if we use compact LDG 
          if( addCompactDG )
          {
            if( problem_.hasCoefficient() )
            {
              CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
              addLiftingOperator(coeffCaller,en,
                                 geo,volQuad,
                                 numDofs,r_e_,matrixEn);

            }
            else 
            {
              CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
              addLiftingOperator(coeffCaller,en,
                                 geo,volQuad,
                                 numDofs,r_e_,matrixEn);
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

    // resize memory for r_e and l_e functions 
    void resizeTemporaryFunctions(const EntityType& en,
                                  TemporaryLocalFunctionArrayType& r_e,
                                  const int numDofs, const int numGradBase) const 
    {
      if( r_e.size() < numDofs )
      {
        r_e.resize(0);
        r_e.resize( numDofs );
        rRets_.resize( numDofs );
        rRetsCoeff_.resize( numDofs );
        for(int i=0; i<numDofs; ++i) 
        {  
          std::auto_ptr< TemporaryLocalFunctionType > 
            ptr(  new TemporaryLocalFunctionType ( gradientSpace_ ) );
          r_e[i] = ptr; 
        }
      }

      if( eta_.size() < numGradBase )
      {
        eta_.resize( numGradBase );
        etaNeigh_.resize( numGradBase );
      }

      for(int i=0; i<numDofs; ++i) 
      {
        TemporaryLocalFunctionType& re = (*r_e[i]);
        re.init ( en );
        assert( re.numDofs() == numGradBase );
        for(int m=0; m<numGradBase; ++m) 
        {
          re[m] = 0;
        }
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

      // get number of base functions for gradient space 
      const int numGradBase = enSet.numBaseFunctions();

      // if we use compact LDG initialize helper functions 
      if( compactLDG_ ) 
      {
        // resize and reset temporary functions 
        resizeTemporaryFunctions( en, r_e_ , numDofs, numGradBase );
#ifndef DG_DOUBLE_FEATURE
        resizeTemporaryFunctions( nb, r_e_neigh_ , numDofs, numGradBase );
#endif
      }

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
        // C_12 stabilization factor 
        const RangeFieldType C_12 = (unitNormal * upwind_ < 0) ? -0.5 : 0.5;
        useInterior = ( C_12 < 0 );
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
          bsetEn.evaluate(k,faceQuadInner[l], phi_[k]);
          // eval gradient for en 
          tau_[k] = bsetEn.evaluateGradientSingle(k, en, faceQuadInner[l] , normEn);  

          // neighbor stuff 
          bsetNeigh.evaluate(k,faceQuadOuter[l], phiNeigh_[k] );      
          // eval gradient for nb 
          tauNeigh_[k] = bsetNeigh.evaluateGradientSingle(k, nb, faceQuadOuter[l] , normNb);      
        }
               
        if(compactLDG_)
        {
          assert( rRets_.size() > 0 );
          GradRangeType& tmp = rRets_[0];

          // get numbre of base functions 
          for(int m=0; m<numGradBase; ++m)
          {  
            // eval base functions 
            enSet.evaluate(m, faceQuadInner[l], tmp );
            // apply unit normal 
            eta_[m] = tmp * unitNormal;

#ifndef DG_DOUBLE_FEATURE
            // neighbor stuff 
            nbSet.evaluate(m, faceQuadOuter[l], tmp ); 
            // apply unit Normal
            etaNeigh_[m] = tmp * unitNormal;
#endif
          }

          for(int k=0; k<numDofs; ++k)
          {
            TemporaryLocalFunctionType& r_e = *(r_e_)[k];
#ifndef DG_DOUBLE_FEATURE
            TemporaryLocalFunctionType& r_e_neigh = *(r_e_neigh_)[k];
#endif
            // calculate [u] 
            // which is the jump of phi 
            RangeType phiDiff (phi_[k]);
            phiDiff -= phiNeigh_[k];
            // scale with integration element 
            phiDiff *= intel;

            // calculate coefficients 
            for(int m=0; m<numGradBase; ++m)
            {
              const RangeFieldType mean = 0.5  * phiDiff * (eta_[m] + etaNeigh_[m]);
              const RangeFieldType jump = C_12 * phiDiff * (eta_[m] - etaNeigh_[m]);
              r_e[m] -= (mean + jump);
#ifndef DG_DOUBLE_FEATURE
              r_e_neigh[m] -= (mean + jump);
#endif
            }
          }
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

              ////////////////////////////////////
              //  C12 stabilization 
              ///////////////////////////////////
              if( compactLDG_)
              {
                typedef typename DiscreteGradientSpaceType :: BaseFunctionSetType BaseFunctionSetType;
                const BaseFunctionSetType enSet = gradientSpace_.baseFunctionSet( en );
                const BaseFunctionSetType nbSet = gradientSpace_.baseFunctionSet( nb );

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
      } // end loop quadrature points 

      if( compactLDG_ )
      {

#ifndef DG_DOUBLE_FEATURE
        if( useInterior ) 
#endif 
        {
          if( problem_.hasCoefficient() )
          {
            CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
            addLiftingOperator(coeffCaller,en,
                               en.geometry(),volQuad,
                               numDofs,r_e_,matrixEn);

          }
          else 
          {
            CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
            addLiftingOperator(coeffCaller,en,
                               en.geometry(),volQuad,
                               numDofs,r_e_,matrixEn);
          }
        }
#ifndef DG_DOUBLE_FEATURE
        else 
        {
          VolumeQuadratureType nbQuad(nb, volumeQuadOrd_);
          if( problem_.hasCoefficient() )
          {
            CoefficientCaller<DiscreteModelCallerType,true,false> coeffCaller; 
            addLiftingOperator(coeffCaller,nb,
                               nb.geometry(),nbQuad,
                               numDofs,r_e_neigh_,matrixNb);

          }
          else 
          {
            CoefficientCaller<DiscreteModelCallerType,false,false> coeffCaller; 
            addLiftingOperator(coeffCaller,nb,
                               nb.geometry(),nbQuad,
                               numDofs,r_e_neigh_,matrixNb);
          }
        }
#endif 
      } // end compactLDG 
    }

    template <class CoeffCallerType>
    void addLiftingOperator(CoeffCallerType& coeffCaller,
                            EntityType& en, 
                            const GeometryType& geo,
                            VolumeQuadratureType& volQuad,
                            const int numDofs, 
                            TemporaryLocalFunctionArrayType& r_e_array,
                            LocalMatrixType& matrixEn) const
    {
      const int volNop = volQuad.nop();
      for (int l = 0; l < volNop ; ++l) 
      {
        // evaluate diffusion coefficient 
        coeffCaller.evaluateCoefficient(caller_, en, volQuad, l, coeffEn_ );

        // calculate integration weight 
        const double intel = volQuad.weight(l)
            * geo.integrationElement(volQuad.point(l));

        for(int k=0; k<numDofs; ++k) 
        {
          // evaluate lifting coefficient function 
          (*r_e_[k]).evaluate(volQuad[l] , rRets_[k] );

          // apply diffusion coefficient 
          coeffCaller.applyCoefficient(coeffEn_, rRets_[k], rRetsCoeff_[k] );
        }

        for(int k=0; k<numDofs; ++k) 
        {
          {
            RangeFieldType val = rRets_[k] * rRetsCoeff_[k];
            val *= intel;
            matrixEn.add(k, k, val);
          }
          for(int j=k+1; j<numDofs; ++j) 
          {
            RangeFieldType val = rRets_[k] * rRetsCoeff_[j];
            val *= intel;

            matrixEn.add(k, j, val);
            matrixEn.add(j, k, val);
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
    GridPartType & gridPart_;
    const DiscreteGradientSpaceType gradientSpace_;

    const LocalIdSetType & localIdSet_;
    const GridWidthType& gridWidth_;
    
    // time
    double time_;

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

    mutable TemporaryLocalFunctionArrayType r_e_;
#ifndef DG_DOUBLE_FEATURE
    mutable TemporaryLocalFunctionArrayType r_e_neigh_;
    DomainType upwind_;
#endif

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
