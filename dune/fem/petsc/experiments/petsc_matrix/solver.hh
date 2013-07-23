#ifndef SOLVER_HH
#define SOLVER_HH

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

class Solver {
public:

    static void solve ( const Mat *matrix, const Vec *rhs, Vec* solution, PetscReal tolerance, 
                        bool multipleProcesses ) {
        if ( multipleProcesses ) {
            return solveParallel( matrix, rhs, solution, tolerance );
        } else {
            return solveSerial( matrix, rhs, solution, tolerance );
        }
    }

    static void solveSerial ( const Mat *matrix, const Vec *rhs, Vec* solution, PetscReal tolerance ) {
        KSP ksp;
        Dune::Petsc::KSPCreate( &ksp );
        Dune::Petsc::KSPSetOperators( ksp, *matrix, *matrix, DIFFERENT_NONZERO_PATTERN );

        PC pc;
        Dune::Petsc::KSPGetPC( ksp, &pc );
        Dune::Petsc::PCSetType( pc, PCJACOBI );
        Dune::Petsc::KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
        Dune::Petsc::KSPSetFromOptions( ksp );

        Dune::Petsc::KSPSolve( ksp, *rhs, *solution );

        Dune::Petsc::KSPDestroy( &ksp );
    }

    static void solveParallel ( const Mat *matrix, const Vec *rhs, Vec* solution, PetscReal tolerance ) {
        KSP ksp;
        Dune::Petsc::KSPCreate( &ksp );
        Dune::Petsc::KSPSetOperators( ksp, *matrix, *matrix, DIFFERENT_NONZERO_PATTERN );

        PC pc;
        Dune::Petsc::KSPGetPC( ksp, &pc );
        Dune::Petsc::PCSetType( pc, PCJACOBI );
        Dune::Petsc::KSPSetTolerances( ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
        Dune::Petsc::KSPSetFromOptions( ksp );

        Dune::Petsc::KSPSolve( ksp, *rhs, *solution );

        Dune::Petsc::KSPDestroy( &ksp );
    }

private:
};

#endif // #if defined HAVE_PETSC

#endif // SOLVER_HH
