#ifndef SOLVER_HH
#define SOLVER_HH

#if defined HAVE_PETSC

#include <dune/fem/petsc/common/petsccommon.hh>

class Solver {
public:

    static void solve ( const Mat *matrix, const Vec *rhs, Vec* solution, PetscReal tolerance ) {
        KSP ksp;
        Dune::Petsc::KSPCreate( &ksp );
        Dune::Petsc::KSPSetOperators( ksp, *matrix, *matrix, DIFFERENT_NONZERO_PATTERN );

        PC pc;
        Dune::Petsc::KSPGetPC( ksp, &pc );
        Dune::Petsc::PCSetType( pc, PCJACOBI );
        Dune::Petsc::KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
        Dune::Petsc::KSPSetFromOptions( ksp );

        PetscInt localSzRhs;
        PetscInt localSzSolution;
        Dune::Petsc::VecGetLocalSize( *rhs, &localSzRhs );
        Dune::Petsc::VecGetLocalSize( *solution, &localSzSolution );
        assert( localSzRhs == localSzSolution );
        //std::cout << "localRhs: " << localSzRhs << ", localSzSolution: " << localSzSolution << std::endl;

        Dune::Petsc::KSPSolve( ksp, *rhs, *solution );

        Dune::Petsc::KSPDestroy( &ksp );
    }

private:
};

#endif // #if defined HAVE_PETSC

#endif // SOLVER_HH
