#ifndef DUNE_FEM_PETSCINVERSEOPERATORS_HH
#define DUNE_FEM_PETSCINVERSEOPERATORS_HH

#include <limits>

#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>

#if HAVE_PETSC
#include <dune/fem/operator/linear/petscoperator.hh>
#include <dune/fem/misc/petsc/petsccommon.hh>
#include <dune/fem/function/petscdiscretefunction.hh>
#include <dune/fem/solver/parameter.hh>

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
    template< class DF, class Op = Dune::Fem::Operator< DF, DF > >
    class PetscInverseOperator;

    template <class DF, class Op >
    struct PetscInverseOperatorTraits
    {
    private:
      typedef typename DF :: DiscreteFunctionSpaceType SpaceType ;
    public:
      typedef DF                                   DiscreteFunctionType;
      typedef Op                                   OperatorType;
      typedef OperatorType                         PreconditionerType;
      typedef PetscDiscreteFunction< SpaceType >   SolverDiscreteFunctionType;
      typedef PetscLinearOperator< DF, DF >        AssembledOperatorType;
      typedef PetscInverseOperator< DF, Op >       InverseOperatorType;
      typedef PetscSolverParameter                 SolverParameterType;
    };


    /** \brief PETSc KSP solver context for PETSc Mat and PETSc Vec */
    template< class DF, class Op >
    class PetscInverseOperator
    : public InverseOperatorInterface< PetscInverseOperatorTraits< DF, Op > >
    {
    protected:
      // monitor function for PETSc solvers
      static PetscErrorCode
      monitor (KSP ksp, PetscInt it, PetscReal rnorm, void *mctx)
      {
        if( Parameter :: verbose ( Parameter::solverStatistics ) )
        {
          std::cout << "PETSc::KSP:  it = "
                    << std::setw(3) << std::left << it
                    << "   res = " << rnorm << std::endl;
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

      typedef PetscInverseOperatorTraits< DF, Op > Traits;
      typedef InverseOperatorInterface< Traits > BaseType;
      friend class InverseOperatorInterface< Traits >;

      enum class PetscSolver {
          cg        = SolverParameter::cg,
          bicgstab  = SolverParameter::bicgstab,
          gmres     = SolverParameter::gmres,
          minres    = SolverParameter::minres,
          bicg      = SolverParameter::bicg,
          preonly   = SolverParameter::preonly,
          kspoptions  = 0
        };


    public:
      static std::vector< int > supportedSolverMethods()
      {
        return std::vector< int > ({
                                     SolverParameter::gmres, // default solver
                                     SolverParameter::cg,
                                     SolverParameter::bicgstab,
                                     SolverParameter::minres,
                                     SolverParameter::bicg,
                                     SolverParameter::preonly });
      }

      static std::vector< int > supportedPreconditionMethods()
      {
        return std::vector< int > ({
                                    SolverParameter::none,   // no preconditioning
                                    SolverParameter::oas,    // Overlapping Additive Schwarz
                                    SolverParameter::gauss_seidel, // SOR with omega = 1
                                    SolverParameter::sor,    // SOR
                                    SolverParameter::ssor,   // symmetric SOR
                                    SolverParameter::jacobi, // Jacobi preconditioning
                                    SolverParameter::ilu,    // ILU preconditioning
                                    SolverParameter::icc     // Incomplete Cholesky factorization
                                  });
      }

      static std::vector<std::string> extraPreconditionMethods()
      {
        return std::vector< std::string > (
                 {"kspoptions", // =  0,   // use command line options -ksp...
                  "hypre",      // = -1,   // Hypre preconditioning
                  "ml",         // = -2,   // ML preconditioner (from Trilinos)
                  "lu",         // = -3,   // LU factorization
                  "pcgamg",     // = -4    // Petsc internal AMG
                 });
      }

      /** \brief This solver does not offer setting preconditioning from outside
       *  \note This needs the implementation of a PCSHELL object to wrap the preconditioner.
       */
      static const bool preconditioningAvailable = false;

      typedef typename BaseType :: SolverDiscreteFunctionType    PetscDiscreteFunctionType;
      typedef typename BaseType :: OperatorType                  OperatorType;
      typedef typename BaseType :: PreconditionerType            PreconditionerType;

      PetscInverseOperator ( const PetscSolverParameter &parameter = PetscSolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {}

      PetscInverseOperator (  const OperatorType &op, const PetscSolverParameter &parameter = PetscSolverParameter(Parameter::container()) )
      : BaseType( parameter )
      {
        bind( op );
      }

      void bind ( const OperatorType &op )
      {
        BaseType :: bind( op );
        initialize( *parameter_ );
      }

      void unbind ()
      {
        BaseType :: unbind();
        ksp_.reset();
      }

      void printTexInfo(std::ostream& out) const
      {
        out << "Solver: " << solverName_ << " eps = " << parameter_->tolerance() ;
        out  << "\\\\ \n";
      }

    protected:
      void initialize ( const PetscSolverParameter& parameter )
      {
        if( !assembledOperator_ )
          DUNE_THROW(NotImplemented,"Petsc solver with matrix free implementations not yet supported!");

        // Create linear solver context
        ksp_.reset( new KSP() );
        const auto& comm = assembledOperator_->domainSpace().gridPart().comm();

        ::Dune::Petsc::KSPCreate( comm, &ksp() );

        // attach Matrix to linear solver context
        Mat& A = assembledOperator_->exportMatrix();
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A, SAME_PRECONDITIONER);
#else
        ::Dune::Petsc::KSPSetOperators( ksp(), A, A );
#endif

        // allow for non-zero initial guess
        ::Dune::Petsc::KSPSetInitialGuessNonzero( ksp(), PETSC_TRUE );

        // set prescribed tolerances
        PetscInt  maxits    = parameter_->maxIterations();
        PetscReal tolerance = parameter_->tolerance();
        PetscReal omega     = parameter_->relaxation();
        if (parameter_->errorMeasure() == 0)
          ::Dune::Petsc::KSPSetTolerances(ksp(), 1e-50, tolerance, PETSC_DEFAULT, maxits);
        else
          ::Dune::Petsc::KSPSetTolerances(ksp(), tolerance, 1e-50, PETSC_DEFAULT, maxits);

        // if special petsc solver parameter exists use that one, otherwise
        // use solverMethod from SolverParameter
        const auto& reader = parameter.parameter();
        PetscSolver kspType = PetscSolver::gmres;
        if( reader.exists("petsc.kspsolver.method") )
        {
          // see PETSc docu for more types
          const std::string kspNames[] = { "default", "cg", "bicgstab", "gmres", "minres", "gradient", "loop", "superlu", "bicg", "preonly"  };
          kspType = static_cast< PetscSolver >( reader.getEnum("petsc.kspsolver.method", kspNames, int(PetscSolver::gmres) ) );
          std::cout << "WARNING: using deprecated parameter 'petsc.kpsolver.method' use "
                    << parameter.keyPrefix() << "method instead\n";
        }
        else
          kspType = static_cast< PetscSolver >(
              parameter.solverMethod (
                supportedSolverMethods(), { "kspoptions" } )
            );

        if (kspType > PetscSolver::kspoptions)
          solverName_ = SolverParameter::solverMethodTable( static_cast< int >( kspType ) );
        else
          solverName_ = "kspoptions";

        //  select linear solver
        switch( kspType )
        {
          case PetscSolver::cg:
            ::Dune::Petsc::KSPSetType( ksp(), KSPCG );
            break;
          case PetscSolver::bicgstab:
            ::Dune::Petsc::KSPSetType( ksp(), KSPBCGS );
            break;
          case PetscSolver::gmres:
            {
              ::Dune::Petsc::KSPSetType( ksp(), KSPGMRES );
              PetscInt restart = 10;
              if( reader.exists("petsc.gmresrestart") )
              {
                restart = reader.getValue<int>("petsc.gmresrestart", restart );
                std::cout << "WARNING: using deprecated parameter 'petsc.gmresrestart' use "
                    << parameter.keyPrefix() << "gmres.restart instead\n";
              }
              else
                restart = parameter.gmresRestart() ;

              ::Dune::Petsc::KSPGMRESSetRestart( ksp(), restart );
              break;
            }
          case PetscSolver::minres:
            ::Dune::Petsc::KSPSetType( ksp(), KSPMINRES );
            break;
          case PetscSolver::bicg:
            ::Dune::Petsc::KSPSetType( ksp(), KSPBICG );
              break;
          case PetscSolver::preonly:
            ::Dune::Petsc::KSPSetType( ksp(), KSPPREONLY );
              break;
          case PetscSolver::kspoptions:
            // setup solver context from database/cmdline options
            ::Dune::Petsc::KSPSetFromOptions( ksp() );
            ::Dune::Petsc::KSPSetUp( ksp() );
            break;
          default:
            DUNE_THROW(InvalidStateException,"PetscInverseOperator: invalid solver choosen." );
        }

        /////////////////////////////////////////////
        //  preconditioning
        /////////////////////////////////////////////

        int pcType = SolverParameter::none;
        if( reader.exists("petsc.preconditioning.method") )
        {
          const std::string pcNames[] = { "default", "none", "asm", "sor", "jacobi", "ilu", "icc", "superlu", "hypre", "ml", "lu" };
          pcType = reader.getEnum("petsc.preconditioning.method", pcNames, 0 );
          std::cout << "WARNING: using deprecated parameter 'petsc.preconditioning.method' use "
                    << parameter.keyPrefix() << ".preconditioning.method instead\n";
          if (pcType >= 8)
            pcType = 7-pcType;  // hypre=-1, ml=-2, lu=-3
        }
        else
        {
          pcType = parameter.preconditionMethod(
                 supportedPreconditionMethods(),
                 extraPreconditionMethods() );
        }

        // setup preconditioning context
        PC pc;
        ::Dune::Petsc::KSPGetPC( ksp(), &pc );

        switch( pcType )
        {
          case 0:
            // don't setup the pc context twice
            if ( kspType != PetscSolver::kspoptions )
            {
              // setup pc context from database/cmdline options
              ::Dune::Petsc::PCSetFromOptions( pc );
              ::Dune::Petsc::PCSetUp( pc );
            }
            break;
          case SolverParameter::none:
            ::Dune::Petsc::PCSetType( pc, PCNONE );
            break;
          case SolverParameter::oas:
            {
              ::Dune::Petsc::PCSetType( pc, PCASM );
              ::Dune::Petsc::PCSetUp( pc );
              break;
            }
          case SolverParameter::gauss_seidel:
            ::Dune::Petsc::PCSetType( pc, PCSOR );
            ::Dune::Petsc::PCSORSetOmega( pc, 1.0 );
            break;
          case SolverParameter::sor:
            ::Dune::Petsc::PCSetType( pc, PCSOR );
            ::Dune::Petsc::PCSORSetOmega( pc, omega );
            break;
          case SolverParameter::ssor:
            ::Dune::Petsc::PCSetType( pc, PCSOR );
            // set symmetric version
            ::Dune::Petsc::PCSORSetSymmetric( pc, SOR_LOCAL_SYMMETRIC_SWEEP );
            ::Dune::Petsc::PCSORSetOmega( pc, omega );
            break;
          case SolverParameter::jacobi:
            ::Dune::Petsc::PCSetType( pc, PCJACOBI );
            break;
          case -1: // PetscPrec::hypre:
            {
#ifdef PETSC_HAVE_HYPRE
              // set with parameter ...petsc.preconditioning.hypre.method
              int hypreType = parameter.hypreMethod();
              std::string hypre;
              if ( hypreType == PetscSolverParameter::boomeramg )
                hypre = "boomeramg";
              else if ( hypreType == PetscSolverParameter::parasails )
                hypre = "parasails";
              else if ( hypreType == PetscSolverParameter::pilut )
                hypre = "pilut";
              else
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid hypre preconditioner choosen." );

              ::Dune::Petsc::PCSetType( pc, PCHYPRE );
              ::Dune::Petsc::PCHYPRESetType( pc, hypre.c_str() );
              ::Dune::Petsc::PCSetUp( pc );
#else // PETSC_HAVE_HYPRE
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with hypre support." );
#endif // PETSC_HAVE_HYPRE
              break;
            }
          case -2: // PetscPrec::ml:
#ifdef PETSC_HAVE_ML
            ::Dune::Petsc::PCSetType( pc, PCML );
#else // PETSC_HAVE_ML
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with ml support." );
#endif // PETSC_HAVE_ML
            break;
          case SolverParameter::ilu:
            {
              if ( MPIManager::size() > 1 )
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: ilu preconditioner does not work in parallel." );

              // get fill-in level
              PetscInt pcLevel;
              if( reader.exists("petsc.preconditioning.levels") )
              {
                pcLevel = reader.getValue<int>("petsc.preconditioning.levels", 0 );
                std::cout << "WARNING: using deprecated parameter 'petsc.preconditioning.levels' use "
                    << parameter.keyPrefix() << "preconditioning.level instead\n";
              }
              else
                pcLevel = parameter.preconditionerLevel() ;

              ::Dune::Petsc::PCSetType( pc, PCILU );
              ::Dune::Petsc::PCFactorSetLevels( pc, pcLevel );
              break;
            }
            ::Dune::Petsc::PCSetType( pc, PCML );
            break;
          case SolverParameter::icc:
            {
#ifdef PETSC_HAVE_ICC
              if ( MPIManager::size() > 1 )
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: icc preconditioner does not worl in parallel." );

              // get fill-in level
              PetscInt pcLevel;
              if( reader.exists("petsc.preconditioning.levels") )
              {
                pcLevel = reader.getValue<int>("petsc.preconditioning.levels", 0 );
                std::cout << "WARNING: using deprecated parameter 'petsc.preconditioning.levels' use "
                    << parameter.keyPrefix() << "preconditioning.level instead\n";
              }
              else
                pcLevel = parameter.preconditionerLevel() ;


              ::Dune::Petsc::PCSetType( pc, PCICC );
              ::Dune::Petsc::PCFactorSetLevels( pc, pcLevel );
#else // PETSC_HAVE_ICC
              DUNE_THROW( InvalidStateException, "PetscInverseOperator: petsc not build with icc support." );
#endif // PETSC_HAVE_ICC
              break;
            }
          case -3: // PetscPrec::lu:
          case SolverParameter::superlu:
            {
              enum class Factorization { petsc = 0, superlu = 1, mumps = 2 };
              Factorization factorType = Factorization::superlu;
              if (pcType != SolverParameter::superlu)
                factorType = static_cast<Factorization>(parameter.superluMethod());

              ::Dune::Petsc::PCSetType( pc, PCLU );

              if ( factorType == Factorization::petsc )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERPETSC );
              else if ( factorType == Factorization::superlu )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERSUPERLU_DIST );
              else if ( factorType == Factorization::mumps )
                ::Dune::Petsc::PCFactorSetMatSolverPackage( pc, MATSOLVERMUMPS );
              else
                DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid factorization package choosen." );

              ::Dune::Petsc::PCSetUp( pc );
              break;
            }
          case -4: // PetscPrec::pcgamg:
            // requires MATRIX_AIJ, i.e. not blocking of entries
            if( parameter.blockedMode() )
            {
              DUNE_THROW(NotImplemented,"PetscInverseOperator: 'pcgamg' requires 'aij' matrix. Set 'petsc.blockedmode' to false!");
            }
            ::Dune::Petsc::PCSetType( pc, PCGAMG );
            break;

          default:
            DUNE_THROW( InvalidStateException, "PetscInverseOperator: invalid preconditioner choosen." );
        }

        // set monitor in verbose mode for all cores
        // (and then check Parameter::verbose locally inside monitor)
        if( parameter.verbose() && Parameter::verbose( Parameter::solverStatistics ) )
        {
          // only print information about solver type and pc type in extended mode
          if( Parameter::verbose( Parameter::extendedStatistics ) )
            ::Dune::Petsc::KSPView( comm, ksp() );

          ::Dune::Petsc::KSPMonitorSet( ksp(), &monitor, PETSC_NULL, PETSC_NULL);
        }
      }

      int apply( const PetscDiscreteFunctionType& arg, PetscDiscreteFunctionType& dest ) const
      {
        // need to have a 'distributed' destination vector for continuous spaces
        if( dest.space().continuous() )
          dest.dofVector().clearGhost();

        // call PETSc solvers
        ::Dune::Petsc::KSPSolve( *ksp_, *arg.petscVec() , *dest.petscVec() );

        // a continuous solution is 'distributed' so need a communication here
        if( dest.space().continuous() )
        {
          dest.communicate();
        }

        // get number of iterations
        PetscInt its ;
        ::Dune::Petsc::KSPGetIterationNumber( *ksp_, &its );
        KSPConvergedReason reason;
        ::Dune::Petsc::KSPGetConvergedReason( *ksp_, &reason );

        bool converged = int(reason) >= 0 ;

        if( parameter_->verbose() && Parameter::verbose( 1 ) )
        {
          // list of reasons:
          // https://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPConvergedReason.html
          if( converged )
            std::cout << "Converged    reason: ";
          else
            std::cout << "**Diverged** reason: ";
          std::cout << reason << "  linear iterations: " << its << std::endl;
        }

        return (converged) ? its : -its;
      }

    protected:
      KSP & ksp () { assert( ksp_ ); return *ksp_; }

      using BaseType :: assembledOperator_;
      using BaseType :: parameter_;
      using BaseType :: iterations_;

      std::unique_ptr< KSP, KSPDeleter > ksp_;   // PETSc Krylov Space solver context

      std::string solverName_;
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // #ifndef DUNE_FEM_PETSCINVERSEOPERATORS_HH
