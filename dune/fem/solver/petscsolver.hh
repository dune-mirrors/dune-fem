#ifndef DUNE_FEM_PETSCSOLVER_HH
#define DUNE_FEM_PETSCSOLVER_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>

#if HAVE_PETSC
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction.hh>

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
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
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

      // destroy solver context
      struct KSPDeleter
      {
        void operator() ( KSP* p ) const
        {
          if( !p )
            return;

          ::Dune::Petsc::KSPDestroy( p );
          delete p;
        }
      };

      // destroy PC context
      struct PCDeleter
      {
        void operator() ( PC *p ) const
        {
          if( !p )
            return;

          ::Dune::Petsc::PCDestroy( p );
          delete p;
        }
      };

    public:
      typedef DF DiscreteFunctionType;
      typedef DiscreteFunctionType  DestinationType;

      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef PetscDiscreteFunction< DiscreteFunctionSpaceType > PetscDiscreteFunctionType;

      typedef PetscLinearOperator< DiscreteFunctionType, DiscreteFunctionType >  AssembledOperatorType;

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
                             double reduction,
                             double absLimit,
                             int maxIter,
                             bool verbose,
                             const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        parameter_( parameter )
      {
        bind( op );
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
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {
        bind( op );
      }

      PetscInverseOperator ( const OperatorType &op,
                             double reduction,
                             double absLimit,
                             const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max()),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {
        bind( op );
      }
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

      PetscInverseOperator ( double reduction, double absLimit, int maxIter, bool verbose, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_( absLimit ),
        maxIter_( maxIter ),
        verbose_( verbose ),
        parameter_( parameter )
      {}

      /** \brief constructor
       *
       *  \param[in] op        mapping describing operator to invert
       *  \param[in] reduction    reduction epsilon
       *  \param[in] absLimit  absolut limit of residual (not used here)
       *  \param[in] maxIter   maximal iteration steps
       */
      PetscInverseOperator ( double reduction, double absLimit, int maxIter, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( maxIter ),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {}

      PetscInverseOperator ( double reduction, double absLimit, const ParameterReader &parameter = Parameter::container() )
      : reduction_( reduction ),
        absLimit_ ( absLimit ),
        maxIter_( std::numeric_limits< int >::max()),
        verbose_( parameter.getValue< bool >( "fem.solver.verbose", false ) ),
        parameter_( parameter )
      {}

      void bind ( const OperatorType &op )
      {
        op_ = &op;
        matrixOp_ = dynamic_cast<const AssembledOperatorType*> ( &op );
        initialize( parameter_ );
      }

      void unbind ()
      {
        op_ = nullptr; matrixOp_ = nullptr;
        ksp_.reset(nullptr);
        pc_.reset(nullptr);
      }

      void initialize ( const ParameterReader &parameter )
      {
        ksp_.reset( new KSP() );

        // Create linear solver context
        ::Dune::Petsc::KSPCreate( &ksp() );

        enum PetscSolver { petsc_cg         = 0,
                           petsc_bicg       = 1,
                           petsc_bicgstab   = 2,
                           petsc_gmres      = 3
                         };

        // see PETSc docu for more types
        const std::string solverNames [] = { "cg" , "bicg", "bicgstab", "gmres" };
        PetscSolver solverType = (PetscSolver) parameter.getEnum("petsc.kspsolver.method", solverNames, 0 );

        // store solver name
        solverName_ = solverNames[ solverType ];

        //  select linear solver
        switch( solverType )
        {
          case petsc_cg:
            {
              ::Dune::Petsc::KSPSetType( ksp(), KSPCG );
              break;
            }
          case petsc_bicg:
            {
              ::Dune::Petsc::KSPSetType( ksp(), KSPBICG ); // this is the BiCG version of PETSc.
              break;
            }
          case petsc_bicgstab:
            {
              ::Dune::Petsc::KSPSetType( ksp(), KSPBCGS ); // this is the BiCG-stab version of PETSc.
              break;
            }
          case petsc_gmres:
            {
              PetscInt restart = parameter.getValue<int>("petsc.gmresrestart", 10 );
              ::Dune::Petsc::KSPSetType( ksp(), KSPGMRES );
              ::Dune::Petsc::KSPGMRESSetRestart( ksp(), restart );
              break;
            }
          default:
            DUNE_THROW(InvalidStateException,"PetscInverseOperator: wrong solver choosen" );
        }


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
        const PetscPcType serialEnd   = petsc_icc;

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
          type = PCHYPRE;
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

        // set preconditioning context
        if( pcType != petsc_none )
        {
          pc_.reset( new PC() );

          // create preconditioning context
          ::Dune::Petsc::PCCreate( &pc() );
          ::Dune::Petsc::PCSetType( pc(), type );
          ::Dune::Petsc::PCFactorSetLevels( pc(), pcLevel );
        }

        if ( pcType == petsc_hypre )
        {
          // set type of HYPRE preconditioner to boomer-amg
          // there are also other preconditioners in this package
          ::Dune::Petsc::PCHYPRESetType( pc(), "boomeramg" );
        }

        // set preconditioning context, if type != none (otherwise problem when solving)
        if( pcType != petsc_none )
          ::Dune::Petsc::KSPSetPC( ksp(), pc() );

        if( pcType == petsc_mumps )
          ::Dune::Petsc::PCFactorSetMatSolverPackage(pc(),MATSOLVERMUMPS);
        if( pcType == petsc_superlu )
          ::Dune::Petsc::PCFactorSetMatSolverPackage(pc(),MATSOLVERSUPERLU_DIST);

        // check if operator is an assembled operator, otherwise we cannot proceed
        if( ! matrixOp_ )
          DUNE_THROW(NotImplemented,"Petsc solver with matrix free implementations not yet supported!");

        // get matrix from linear operator
        Mat& A = const_cast<Mat &> (matrixOp_->petscMatrix());

        // set operator to PETSc solver context
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A, SAME_PRECONDITIONER);
#else
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A );
#endif
        // set prescribed tolerances
        PetscInt  maxits = maxIter_ ;
        PetscReal reduc  = reduction_;
        ::Dune::Petsc::KSPSetTolerances(ksp(), reduc, 1.e-50, PETSC_DEFAULT, maxits);

        // set monitor in verbose mode
        if( verbose_ )
        {
          ::Dune::Petsc::KSPView( ksp() );
          ::Dune::Petsc::KSPMonitorSet( ksp(), &monitor, PETSC_NULL, PETSC_NULL);
        }
      }

      void setMaxIterations ( std::size_t maxIter )
      {
        maxIter_ = maxIter;
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
      template <class DiscreteFunction>
      void apply( const DiscreteFunction& arg, DiscreteFunction& dest ) const
      {
        // copy discrete functions
        PetscDiscreteFunctionType Arg("PetscSolver::arg", arg.space() );
        Arg.assign( arg );

        // also copy initial destination in case this is used a solver init value
        PetscDiscreteFunctionType Dest("PetscSolver::dest", dest.space() );
        Dest.assign( dest );

        apply( Arg, Dest );
        // copy destination back
        dest.assign( Dest );
      }

      void apply( const PetscDiscreteFunctionType& arg, PetscDiscreteFunctionType& dest ) const
      {
        // call PETSc solvers
        ::Dune::Petsc::KSPSolve( *ksp_, *arg.petscVec() , *dest.petscVec() );

        // for continuous solution we need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        // get number of iterations
        PetscInt its ;
        ::Dune::Petsc::KSPGetIterationNumber( *ksp_, &its );
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

      template <class DiscreteFunction>
      void operator() ( const DiscreteFunction& arg, DiscreteFunction& dest ) const
      {
        apply(arg,dest);
      }

    protected:
      PC & pc () { assert( pc_ ); return *pc_; }
      KSP & ksp () { assert( ksp_ ); return *ksp_; }

      const OperatorType  *op_ = nullptr; // linear operator
      const AssembledOperatorType* matrixOp_ = nullptr; // assembled operator

      std::unique_ptr< KSP, KSPDeleter > ksp_;   // PETSc Krylov Space solver context
      std::unique_ptr< PC, PCDeleter > pc_;    // PETSc perconditioning context
      double reduction_;
      double absLimit_;
      int maxIter_;
      bool verbose_ ;
      mutable int iterations_ = 0;
      ParameterReader parameter_;
      std::string solverName_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCSOLVER_HH
