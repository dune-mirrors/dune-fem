#ifndef DUNE_FEM_PETSCAVAILABLE_HH
#define DUNE_FEM_PETSCAVAILABLE_HH

#include <limits>

#include <dune/fem/io/parameter.hh>
#include <dune/fem/solver/inverseoperatorinterface.hh>
#include <dune/fem/solver/parameter.hh>

namespace Dune
{

  namespace Fem
  {

    //=====================================================================
    // Implementation of available solver/preconditioners for
    // PETSc matrix based Krylov solvers
    //=====================================================================

    /** @ingroup OEMSolver
        @{
    **/

    /** \brief PETSc KSP solver context for PETSc Mat and PETSc Vec */
    class PetscInverseOperatorAvailable
    {
    protected:
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
    };

  ///@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_PETSCINVERSEOPERATORS_HH
