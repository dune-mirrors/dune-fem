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
                             bool verbose,
                             const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        iterations_( 0 ),
        solverName_("none"),
        precondName_("none")
      {
        initialize( parameter );
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
                             int maxIter,
                             const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        solverName_("none"),
        precondName_("none")
      {
        initialize( parameter );
      }

      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             const ParameterReader &parameter = Parameter::container() )
      : op_( op ),
        reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max()),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        iterations_( 0 ),
        solverName_("none"),
        precondName_("none")
      {
        initialize( parameter );
      }

      //! destructor freeing KSP context
      ~PetscInverseOperator()
      {
        // destroy PC context
        ::Dune::Petsc::PCDestroy( &pc_ );
        // destroy solver context
        ::Dune::Petsc::KSPDestroy( &ksp_ );
      }

      void initialize ( const ParameterReader &parameter )
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
        PetscSolver solverType = (PetscSolver) parameter.getEnum("petsc.kspsolver.method", solverNames, 0 );

        //  select linear solver
        if ( solverType == petsc_cg )
          ::Dune::Petsc::KSPSetType( ksp_, KSPCG );
        else if ( solverType == petsc_bicg )
          ::Dune::Petsc::KSPSetType( ksp_, KSPBICG ); // this is the BiCG version of PETSc.
        else if ( solverType == petsc_bicgstab )
          ::Dune::Petsc::KSPSetType( ksp_, KSPBCGS ); // this is the BiCG-stab version of PETSc.
        else if ( solverType == petsc_gmres )
        {
          PetscInt restart = parameter.getValue<int>("petsc.gmresrestart", 10 );
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

        enum PetscPcType { petsc_none       = 0,   // no preconditioning
                           // parallel preconditioners
                           petsc_asm        = 1,   // Additive Schwarz
                           petsc_sor        = 2,   // SOR and SSOR
                           petsc_jacobi     = 3,   // Jacobi preconditioning
                           // requiring additional packages
                           petsc_hypre      = 4,   // Hypre preconditioning
                           petsc_ml         = 5,   // ML preconditioner (from Trilinos)
                           // serial preconditioners
                           petsc_ilu        = 6,   // ILU preconditioning
                           petsc_lu         = 7,   // LU factorization
                           petsc_icc        = 8,   // Incomplete Cholesky factorization
                           // direct solvers
                           petsc_mumps      = 9,   // use mumps
                           petsc_superlu    = 10,   // use superlu-dist
                           end              = 11
                         };
        const PetscPcType serialStart = petsc_ilu;
        const PetscPcType serialEnd = petsc_icc;

        PCType type = PCNONE ;

        // see PETSc docu for more types
        const std::string preconNames [end] = { "none", "asm", "sor", "jacobi",
                                                "hypre", "ml",
                                                "ilu-n", "lu", "icc",
                                                "mumps", "superlu" };
        PetscPcType pcType = (PetscPcType) parameter.getEnum("petsc.preconditioning.method", preconNames, 0 );

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
          type = PCHYPRE;
        }
        else if ( pcType == petsc_ml )
          type = PCML;
        else if ( pcType == petsc_ilu )
          type = PCILU;
        else if( pcType == petsc_lu )
          type = PCLU;
        else if( pcType == petsc_icc )
          type = PCICC;
        else if( pcType == petsc_mumps )
        {
          /* could use
          KSPSetType(ksp,KSPPREONLY);
          */
          type = PCLU;
        }
        else if( pcType == petsc_superlu )
        {
          /* could use
          KSPSetType(ksp,KSPPREONLY);
          */
          type = PCLU;
        }
        else
          DUNE_THROW(InvalidStateException,"PetscInverseOperator: wrong preconditiong choosen" );

        // check whether the preconditioner can be used in parallel
        if( MPIManager :: size() > 1 && pcType >= serialStart && pcType <= serialEnd )
        {
          if( MPIManager :: rank() == 0 )
          {
            std::cerr << "WARNING: PetscInverseOperator: " << preconNames[ pcType ]  << " preconditioning does not work in parallel, falling back to Additive Schwarz!!" << std::endl;
          }
          type = PCASM;
        }

        // get level of perconditioner iterations
        PetscInt pcLevel = 1 + parameter.getValue<int>("petsc.preconditioning.iterations", 0 );

        // create preconditioning context
        ::Dune::Petsc::PCCreate( &pc_ );
        ::Dune::Petsc::PCSetType( pc_, type );
        ::Dune::Petsc::PCFactorSetLevels( pc_, pcLevel );

        if ( pcType == petsc_hypre )
        {
          // set type of HYPRE preconditioner to boomer-amg
          // there are also other preconditioners in this package
          ::Dune::Petsc::PCHYPRESetType( pc_, "boomeramg" );
        }

        // set preconditioning context
        ::Dune::Petsc::KSPSetPC( ksp_, pc_ );
        if( pcType == petsc_mumps )
          ::Dune::Petsc::PCFactorSetMatSolverPackage(pc_,MATSOLVERMUMPS);
        if( pcType == petsc_superlu )
          ::Dune::Petsc::PCFactorSetMatSolverPackage(pc_,MATSOLVERSUPERLU_DIST);

        // get matrix from linear operator
        Mat& A = const_cast< Mat & > (op_.petscMatrix());

        // set operator to PETSc solver context
        // ::Dune::Petsc::KSPSetOperators( ksp_, A, A, DIFFERENT_NONZERO_PATTERN);
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
        ::Dune::Petsc::KSPSetOperators( ksp_, A, A, SAME_PRECONDITIONER);
#else
        ::Dune::Petsc::KSPSetOperators( ksp_, A, A );
#endif
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
      {}

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

        // for continuous solution we need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

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

      PetscInverseOperator () = delete;
      PetscInverseOperator( const PetscInverseOperator& ) = delete;
      PetscInverseOperator& operator= ( const PetscInverseOperator& ) = delete;

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
