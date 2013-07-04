#ifndef SCHEME_HH
#define SCHEME_HH

#include <assert.h>
#include <iostream>
#include <cmath>

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

#include "function.hh"
#include "l2error.hh"
#include "laplaceassembler.hh"
#include "solver.hh"
#include "printer.hh"

template< typename S >
class Scheme {
public:
    typedef S ScalarType;
    

    Scheme ( unsigned int initialGridSize,
             unsigned int eocSteps,
             const Function< ScalarType > &f,
             const Function< ScalarType > &exactSolution,
             bool multipleProcesses )
    : initialGridSize_( initialGridSize ),
      eocSteps_( eocSteps ),
      f_( f ),
      exactSolution_( exactSolution ),
      multipleProcesses_( multipleProcesses )
    {
        assert( f_.lowerBound() == exactSolution.lowerBound() );
        assert( f_.upperBound() == exactSolution.upperBound() );
    }


    void performAllSteps () const {

        int myRank; 
        MPI_Comm_rank( PETSC_COMM_WORLD, &myRank ); 
    
        unsigned int currentGridSize = initialGridSize_;
        ScalarType newError = 0.;
        ScalarType oldError = performEOCStep( currentGridSize, false );
        ScalarType eoc;
        ScalarType log2 = std::log( 2 );

        for ( unsigned int step = 1; step < eocSteps_; ++step ) {

            currentGridSize *= 2;
            newError = performEOCStep( currentGridSize, false );
            eoc = std::log( oldError/newError ) / log2;

            if( myRank == 0 ) {
                std::cout << "\n\n+++++ In EOC step " << step << " ++++++++++\nOld error: " << oldError 
                          << ", new error: " << newError << "\nEOC: " << eoc << std::endl;
            }

            oldError = newError;
        }
    }

    // performs an eoc step and returns the error
    ScalarType performEOCStep ( unsigned int gridSize, bool printValues = true ) const {

    
        // global settings
        Mat A;
        Vec rhs;
        Vec solution;

        LaplaceAssembler assembler( gridSize, f_.lowerBound(), f_.upperBound(), multipleProcesses_ );
        
        assembler.assembleMatrix( &A );
        assembler.assembleRHS( f_, &rhs );
        assembler.initializeSolutionVector( &solution );

        //PETSCErrorCheck( VecView( rhs, PETSC_VIEWER_STDOUT_WORLD ) );
        //PETSCErrorCheck( MatView( A, PETSC_VIEWER_STDOUT_WORLD ) );

        Solver::solve( &A, &rhs, &solution, static_cast< PetscReal >( 1.e-7 ), multipleProcesses_ );

        if ( printValues ) {
            Printer::printCout( gridSize, solution, exactSolution_, multipleProcesses_ );
        }

        PetscScalar error = L2Error< PetscScalar >::calculateError( solution, exactSolution_, gridSize );

        // finalize the program
        Dune::Petsc::VecDestroy( &rhs );
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
    bool multipleProcesses_;

};

#endif // #if defined HAVE_PETSC

#endif // SCHEME_HH
