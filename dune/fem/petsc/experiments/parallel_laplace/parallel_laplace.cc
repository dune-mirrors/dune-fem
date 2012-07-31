#include <config.h>

#include <iostream>

#include <mpi.h>

#include <dune/fem/petsc/common/petsccommon.hh>

#include "function.hh"
#include "scheme.hh"

int main ( int argc, char **argv ) {
    Dune::Petsc::PetscInitialize( &argc, &argv, static_cast< char* >( 0 ), "A small testing program" );

    // only run in parallel
    int commSize;
    MPI_Comm_size( MPI_COMM_WORLD, &commSize );
    if( commSize <= 1 ) {
        std::cerr << "please run '" << argv[ 0 ] << "' with at least 2 MPI jobs!\n";
        Dune::Petsc::PetscFinalize();
        return 1;
    } 

    PetscScalar a = -1.;
    PetscScalar b = 1.;
    unsigned int initialGridSize = 10;
    FunctionSine< PetscScalar > rhsFunc( a, b );
    SolutionSine< PetscScalar > exactSolution( a, b );

    Scheme< PetscScalar > scheme( initialGridSize, 10, rhsFunc, exactSolution );

    if( argc > 1 && std::string( argv[ 1 ] ) == "eoc" )  {
        scheme.performAllSteps();
    } else {
        std::cerr << "L2 error: " << scheme.performEOCStep( initialGridSize, false ) << std::endl;
    }

    Dune::Petsc::PetscFinalize();
    return 0;
}
