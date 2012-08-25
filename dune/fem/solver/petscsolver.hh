#ifndef DUNE_FEM_PETSCSOLVER_HH 
#define DUNE_FEM_PETSCSOLVER_HH 

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#if HAVE_PETSC 
#include <dune/fem/misc/petsc/petsccommon.hh>

namespace Dune
{

  namespace Fem 
  {

    //=====================================================================
    // Implementation for PETSc matrix based Krylov solvers 
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    // PETScSolver 
    // --------------

    /** \brief PETSc KSP solver context for PETSc Mat and PETSc Vec */
    template< class DF, class Op >
    class PetscInverseOperator 
      : public Operator< DF, DF >
    {
    protected:
      // monitor function for PETSc solvers 
      static PetscErrorCode 
      monitor (KSP ksp, int it, PetscReal rnorm, void *mctx)
      {
        if( Parameter :: verbose () )
        {
          std::cout << "PETSc::KSP:  it = " << it << "    res = " << rnorm << std::endl;
        }
        return PetscErrorCode(0);
      }

    public:  
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;
      typedef Op OperatorType;

    public:
      /** \brief constructor
       *
       *  \param[in] op Mapping describing operator to invert
       *  \param[in] reduction reduction epsilon
       *  \param[in] absLimit absolut limit of residual (not used here)
       *  \param[in] maxIter maximal iteration steps
       *  \param[in] verbose verbosity
       *
       *  \note PETSc Krylov solvers uses the relative reduction. 
       */
      PetscInverseOperator ( const OperatorType &op,
                             double  reduction,
                             double absLimit,
                             int maxIter,
                             bool verbose )
      : op_( op ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        solverName_("none"),
        precondName_("none")
      {
        initialize();
      }

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             int maxIter = std::numeric_limits< int >::max() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( Parameter::getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        solverName_("none"),
        precondName_("none")
      {
        initialize();
      }

      //! destructor freeing KSP context 
      ~PetscInverseOperator() 
      {
        // destroy PC context
        ::Dune::Petsc::PCDestroy( &pc_ );
        // destroy solver context 
        ::Dune::Petsc::KSPDestroy( &ksp_ );
      }

      void initialize () 
      {
        // Create linear solver context
        ::Dune::Petsc::KSPCreate( &ksp_ );
    
        enum PetscSolver { petsc_cg         = 0,
                           petsc_bicg       = 1,
                           petsc_bicgstab   = 2,
                           petsc_gmres      = 3 
                         };

        // see PETSc docu for more types 
        const std::string solverNames [] = { "cg" , "bicg", "bicgstab", "gmres" };
        PetscSolver solverType = (PetscSolver) Parameter :: getEnum("petsc.kspsolver.method", solverNames, 0 );

        //  select linear solver 
        if ( solverType == petsc_cg )
          ::Dune::Petsc::KSPSetType( ksp_, KSPCG );
        else if ( solverType == petsc_bicg )
          ::Dune::Petsc::KSPSetType( ksp_, KSPBICG ); // this is the BiCG version of PETSc.
        else if ( solverType == petsc_bicgstab )
          ::Dune::Petsc::KSPSetType( ksp_, KSPBCGS ); // this is the BiCG-stab version of PETSc.
        else if ( solverType == petsc_gmres )
        {
          PetscInt restart = Parameter :: getValue<int>("petsc.gmresrestart", 10 );
          ::Dune::Petsc::KSPSetType( ksp_, KSPGMRES );
          ::Dune::Petsc::KSPGMRESSetRestart( ksp_, restart ); 
        }
        else
        {
          DUNE_THROW(InvalidStateException,"PetscInverseOperator: wrong solver choosen" );
        }

        // store solver name 
        solverName_ = solverNames[ solverType ];

        /////////////////////////////////////////////
        //  preconditioning 
        /////////////////////////////////////////////

        enum PetscPcType { petsc_none     = 0, // no preconditioning 
                           // parallel preconditioners 
                           petsc_asm      = 1, // Additive Schwarz 
                           petsc_sor      = 2, // SOR and SSOR 
                           petsc_jacobi   = 3, // Jacobi preconditioning 
                           petsc_hypre    = 4, // Hypre preconditioning
                           // serial preconditioners 
                           petsc_ilu      = 5, // ILU preconditioning 
                           petsc_lu       = 6, // LU factorization 
                           petsc_icc      = 7,  // Incomplete Cholesky factorization
                         };

        const PCType type = PCNONE ;

        // see PETSc docu for more types 
        const std::string preconNames [] = { "none", "asm", "sor", "jacobi", "hypre",  "ilu-n", "lu", "icc" };
        PetscPcType pcType = (PetscPcType) Parameter :: getEnum("petsc.preconditioning.method", preconNames, 0 );

        if( pcType == petsc_none )
          type = PCNONE ;
        else if( pcType == petsc_asm ) 
          type = PCASM;
        else if( pcType == petsc_sor ) 
          type = PCSOR;
        else if ( pcType == petsc_jacobi ) 
          type = PCJACOBI;
        else if ( pcType == petsc_hypre ) 
        {
          #if defined HAVE_PETSC_HYPRE
          type = PCHYPRE;
          #else
          DUNE_THROW(InvalidStateException,"PetscInverseOperator: Trying to use HYPRE preconditioning but PETSc has not been compiled with HYPRE support" );
          #endif // defined HAVE_PETSC_HYPRE
        }
        else if ( pcType == petsc_ilu ) 
          type = PCILU; 
        else if( pcType == petsc_lu ) 
          type = PCLU;
        else if( pcType == petsc_icc ) 
          type = PCICC;
        else 
          DUNE_THROW(InvalidStateException,"PetscInverseOperator: wrong preconditiong choosen" );

        // check whether the preconditioner can be used in parallel
        if( MPIManager :: size() > 1 && pcType > petsc_hypre ) 
        {
          if( MPIManager :: rank() == 0 )
          {
            std::cerr << "WARNING: PetscInverseOperator: " << preconNames[ pcType ]  << " preconditioning does not work in parallel, falling back to Additive Schwarz!!" << std::endl; 
          }
          type = PCASM;
        }

        // get level of perconditioner iterations 
        PetscInt pcLevel = 1 + Parameter :: getValue<int>("petsc.preconditioning.iterations", 0 );

        // create preconditioning context 
        ::Dune::Petsc::PCCreate( &pc_ );
        ::Dune::Petsc::PCSetType( pc_, type );
        ::Dune::Petsc::PCFactorSetLevels( pc_, pcLevel );

        // :: PCHYPRESetType( pc_, "boomeramg" );

        // set preconditioning context 
        ::Dune::Petsc::KSPSetPC( ksp_, pc_ );

        // get matrix from linear operator 
        Mat& A = const_cast< Mat & > (op_.petscMatrix());

        // set operator to PETSc solver context 
        ::Dune::Petsc::KSPSetOperators( ksp_, A, A, DIFFERENT_NONZERO_PATTERN);

        // set prescribed tolerances 
        PetscInt  maxits = maxIter_ ;
        PetscReal reduc  = reduction_;
        ::Dune::Petsc::KSPSetTolerances(ksp_, reduc, 1.e-50, PETSC_DEFAULT, maxits);

        // set monitor in verbose mode 
        if( verbose_ ) 
        {
          ::Dune::Petsc::KSPView( ksp_ );
          ::Dune::Petsc::KSPMonitorSet( ksp_, &monitor, PETSC_NULL, PETSC_NULL);
        }
      }

      void prepare (const DiscreteFunctionType& Arg, DiscreteFunctionType& Dest) const
      {
      }

      void finalize () const
      {}

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: " << solverName_ << " eps = " << reduction_ ;
        out  << "\\\\ \n";
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void apply( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        // call PETSc solvers 
        ::Dune::Petsc::KSPSolve(ksp_, *arg.petscVec() , *dest.petscVec() );

        // get number of iterations 
        PetscInt its ;
        ::Dune::Petsc::KSPGetIterationNumber( ksp_, &its );
        iterations_ = its;
      }

      // return number of iterations 
      int iterations() const 
      {
        return iterations_;
      }

      //! return accumulated communication time
      double averageCommTime() const 
      {
        return -1;
      }

      /** \brief solve the system 
          \param[in] arg right hand side 
          \param[out] dest solution 
      */
      void operator() ( const DiscreteFunctionType& arg, DiscreteFunctionType& dest ) const
      {
        apply(arg,dest);
      }

    private:
      // no copying 
      PetscInverseOperator ();
      PetscInverseOperator( const PetscInverseOperator& ) ;
      PetscInverseOperator& operator= ( const PetscInverseOperator& );

    protected:
      const OperatorType &op_; // linear operator 
      KSP ksp_;   // PETSc Krylov Space solver context 
      PC  pc_;    // PETSc perconditioning context 
      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_;
      std::string solverName_;
      std::string precondName_;
    }; 

  ///@}

  } // namespace Fem 

} // namespace Dune 

#endif // #if HAVE_PETSC 

#endif // #ifndef DUNE_FEM_PETSCSOLVER_HH
