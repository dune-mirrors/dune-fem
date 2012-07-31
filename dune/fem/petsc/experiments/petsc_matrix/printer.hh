#ifndef PRINTER_HH
#define PRINTER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>

#if defined HAVE_PETSC

#include <dune/fem/petsc/common/petsccommon.hh>

#include "function.hh"


class Printer {
public:
    
    // prints both solutions to std::cout in a format that gnuplot understands
    template< typename S >
    static void printCout ( unsigned int gridSize,
                            Vec numericalSolution,
                            const Function< S > &exactSolution,
                            bool multipleProcesses = false )
    {
        PetscInt vecSize;
        Dune::Petsc::VecGetSize( numericalSolution, &vecSize );
        assert( static_cast< unsigned int >( vecSize ) == gridSize );
        
        PetscScalar a = exactSolution.lowerBound();
        PetscScalar b = exactSolution.upperBound();
        PetscScalar h = (b-a) / static_cast< PetscScalar >( gridSize - 1 );
        //std::cout << "h in Printer: " << h << std::endl;
        //std::cout << "a in Printer: " << a << std::endl;
        //std::cout << "b in Printer: " << b << std::endl;

        /*
         * VecGetValues currently only supports coordinate access for the coordinates of 
         * this process. We use a custom viewer and print all the coordinates to a temporary file...
         */
        const char* bufferFileName = "tmpVecBuffer611e08b4f16a5041f339d76d0ea101da.out";
        PetscViewer viewer;
        MPI_Comm comm;
        if ( multipleProcesses ) {
            comm = PETSC_COMM_WORLD;
        } else {
            comm = PETSC_COMM_SELF;
        }
        Dune::Petsc::PetscViewerASCIIOpen( comm, bufferFileName, &viewer );
        Dune::Petsc::VecView( numericalSolution, viewer );

        std::ifstream vecFile( bufferFileName );
        std::string vecFileLine;

        for( unsigned int i=0; i < gridSize; ++i ) {
            PetscScalar currentX = a + i*h;
            vecFile >> vecFileLine;
            std::cout << currentX << "\t" << vecFileLine << "\t\t" << exactSolution.evaluate( currentX ) << "\n";
        }

        // clean up, i.e. remove the temporary file
        std::remove( bufferFileName );

    }

private:
    // prohibited methods
    Printer ();

    #if 0
    PetscScalar a_;
    PetscScalar b_;
    unsigned int gridSize_;
    #endif
};

#endif // #if defined HAVE_PETSC

#endif // PRINTER_HH
