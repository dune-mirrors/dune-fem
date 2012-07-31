#include <config.h>
// sys includes
#include <iostream>
#include <string>

// PETSc includes
#include "petscvec.h"
#include "petscmat.h"
#include "petscsys.h"

// my includes
#include "function.hh"
#include "laplaceassembler.hh"
#include "solver.hh"
#include "printer.hh"
#include "l2error.hh"
#include "scheme.hh"

/*
 * This shall be a FEM-scheme for 1D...
 */

bool multipleProcesses = true;

int main ( int argc, char **argv ) {
    PetscInitialize( &argc, &argv, static_cast< char* >( 0 ), "A small testing program" );

    PetscScalar a = -1.;
    PetscScalar b = 1.;
    unsigned int initialGridSize = 100;
    FunctionSine< PetscScalar > rhsFunc( a, b );
    SolutionSine< PetscScalar > exactSolution( a, b );

    Scheme< PetscScalar > scheme( initialGridSize, 10, rhsFunc, exactSolution, multipleProcesses );

    if ( argc > 1 && std::string( argv[ 1 ] ) == "eoc" )  {
        scheme.performAllSteps();
    } else {
        std::cerr << "L2 error: " << scheme.performEOCStep( initialGridSize, false ) << std::endl;
    }

    PetscFinalize();
    return 0;
}
