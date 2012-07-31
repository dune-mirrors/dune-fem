#ifndef SCHEME_HH
#define SCHEME_HH

#include <assert.h>
#include <iostream>
#include <sstream>
#include <cmath>

#if defined HAVE_PETSC

#include <dune/fem/petsc/common/petsccommon.hh>

#include "petscmat.h"
#include "petscvec.h"

#include "solver.hh"
#include "printer.hh"
#include "function.hh"
#include "l2error.hh"
#include "laplaceassembler.hh"

template< typename S >
class Scheme {
public:
    typedef S ScalarType;
    

    Scheme ( unsigned int initialGridSize,
             unsigned int eocSteps,
             const Function< ScalarType > &f,
             const Function< ScalarType > &exactSolution )
    : initialGridSize_( initialGridSize ),
      eocSteps_( eocSteps ),
      f_( f ),
      exactSolution_( exactSolution )
    {
        assert( f_.lowerBound() == exactSolution.lowerBound() );
        assert( f_.upperBound() == exactSolution.upperBound() );
    }


    void performAllSteps () const {

        int myRank; 
        MPI_Comm_rank( PETSC_COMM_WORLD, &myRank ); 
    
        unsigned int currentGridSize = initialGridSize_;
        ScalarType newError = 0.;
        ScalarType oldError = performEOCStep( currentGridSize, false, 0 );
        ScalarType eoc;
        ScalarType log2 = std::log( 2 );

        for( unsigned int step = 1; step < eocSteps_; ++step ) {

            currentGridSize *= 2;
            newError = performEOCStep( currentGridSize, false, step );
            eoc = std::log( oldError/newError ) / log2;

            if( myRank == 0 ) {
                std::cout << "\n\n+++++ In EOC step " << step << " ++++++++++\nOld error: " << oldError 
                          << ", new error: " << newError << "\nEOC: " << eoc
                          << ", current grid size: " << currentGridSize << std::endl;
            }

            oldError = newError;
        }
    }

    // performs an eoc step and returns the error
    ScalarType performEOCStep ( unsigned int gridSize, bool printValues = true, int eocStep = 0 ) const {

    
        // global settings
        Mat A;
        Vec rhs;
        Vec solution;

        LaplaceAssembler assembler( gridSize, f_.lowerBound(), f_.upperBound() );
        
        assembler.assembleRHS( f_, &rhs, eocStep );
        assembler.assembleMatrix( &A );
        assembler.initializeSolutionVector( &solution );

        //PETSCErrorCheck( VecView( rhs, PETSC_VIEWER_STDOUT_WORLD ) );
        //PETSCErrorCheck( MatView( A, PETSC_VIEWER_STDOUT_WORLD ) );

        Solver::solve( &A, &rhs, &solution, static_cast< PetscReal >( 1.e-7 ) );

        
        std::stringstream gnuplotFilename;
        gnuplotFilename << "solution.out.gnuplot." << eocStep;
        std::ofstream outFile( gnuplotFilename.str() );
        Printer::printVec( gridSize, f_.lowerBound(), f_.upperBound(), solution, outFile );


        PetscScalar error = L2Error< PetscScalar >::calculateError( solution, exactSolution_, gridSize, eocStep );

        // finalize the program
        assembler.finalizeRHS( &rhs );
        Dune::Petsc::VecDestroy( &solution );
        Dune::Petsc::MatDestroy( &A );

        return static_cast< ScalarType >( error );
    }
    


private:
    // prohibited methods
    Scheme ();

    // data fields
    unsigned int initialGridSize_;
    unsigned int eocSteps_;
    const Function< S > &f_;
    const Function< S > &exactSolution_;

};

#endif // #if defined HAVE_PETSC

#endif // SCHEME_HH
