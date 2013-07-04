#ifndef LAPLACEASSEMBLER_HH
#define LAPLACEASSEMBLER_HH

#include <assert.h>
#include <iostream>

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

#include "function.hh"

/*
 * Class for PETSc matrix and vector assembly
 */
//template< bool multipleProcesses >
class LaplaceAssembler {
public:

    LaplaceAssembler ( unsigned int gridSize, PetscScalar a, PetscScalar b, bool multipleProcesses )
    : gridSize_( gridSize ),
      a_( a ),
      b_( b ),
      h_( static_cast< PetscScalar >( b_ - a_ )/ static_cast< PetscScalar >( gridSize_ - 1 ) ),
      hSquare_( h_ * h_ ),
      multipleProcesses_( multipleProcesses )
    {
        //std::cout << "h in Assembler: " << h_ << std::endl;
        assert( b_ > a_ );
    }

    // Creates and assembles the matrix 'matrix'. Calls, among other functions, MatCreate(...) on
    // 'matrix', so 'matrix' should be supplied entirely uninitialized...
    void assembleMatrix ( Mat *matrix ) const {
        if ( multipleProcesses_ ) {
            return assembleMatrixParallel( matrix );
        } else {
            return assembleMatrixSerial( matrix );
        }
    }

    // assemble the right hand side vector
    template< typename S >
    void assembleRHS ( const Function< S > &f, Vec *vector ) const {
        if ( multipleProcesses_ ) {
            return assembleRHSParallel( f, vector );
        } else {
            return assembleRHSSertial( f, vector );
        }
    }

    void initializeSolutionVector ( Vec *vector ) const {
        if ( multipleProcesses_ ) {
            return initializeSolutionVectorParallel( vector );
        } else {
            return initializeSolutionVectorSerial( vector );
        }
    }

private:
    /*
     * private methods
     */

    void assembleMatrixSerial ( Mat *matrix ) const {

        Dune::Petsc::MatCreate( matrix );
        Dune::Petsc::MatSetSizes( *matrix, PETSC_DECIDE, PETSC_DECIDE, gridSize_, gridSize_ );
        Dune::Petsc::MatSetFromOptions( *matrix );

        /* 
         * assemble the matrix elements
         */
        PetscInt i;
        PetscScalar values[ 3 ]; values[ 0 ] = -1.; values[ 1 ] = 2.; values[ 2 ] = -1.;
        PetscInt colIndices[ 3 ];
        // set all but the first and last rows
        for ( i=1; i < static_cast< PetscInt >( gridSize_ ) - 1; ++i ) {
            colIndices[ 0 ] = i-1; colIndices[ 1 ] = i; colIndices[ 2 ] = i+1;
            Dune::Petsc::MatSetValues( *matrix, 1, &i, 3, colIndices, values, INSERT_VALUES );
        }
        // set the last row
        colIndices[ 0 ] = static_cast< PetscInt >( gridSize_ ) - 2; 
        colIndices[ 1 ] = static_cast< PetscInt >( gridSize_ ) - 1;
        values[ 0 ] = 0.; values[ 1 ] = 1.;
        Dune::Petsc::MatSetValues( *matrix, 1, &i, 2, colIndices, values, INSERT_VALUES );
        // set the first row
        i = 0;
        colIndices[ 0 ] = 0; colIndices[ 1 ] = 1;
        values[ 0 ] = 1.; values[ 1 ] = 0.;
        Dune::Petsc::MatSetValues( *matrix, 1, &i, 2, colIndices, values, INSERT_VALUES );

        // MatAssemblyBegin should be called _after_ all MatSetValues calls...
        Dune::Petsc::MatAssemblyBegin( *matrix, MAT_FINAL_ASSEMBLY );
        Dune::Petsc::MatAssemblyEnd( *matrix, MAT_FINAL_ASSEMBLY );
    }

    void assembleMatrixParallel ( Mat *matrix ) const {

        Dune::Petsc::MatCreate( matrix );
        Dune::Petsc::MatSetSizes( *matrix, PETSC_DECIDE, PETSC_DECIDE, gridSize_, gridSize_ );
        Dune::Petsc::MatSetFromOptions( *matrix );

        // Do I need this? Is this only for speedup?
        /*
        PETSCErrorCheck( MatMPIAIJSetPreallocation( *matrix, 5, PETSC_NULL, 5, PETSC_NULL ) );
        PETSCErrorCheck( MatSeqAIJSetPreallocation( *matrix, 5, PETSC_NULL ) );
        */

        PetscInt start, end;
        Dune::Petsc::MatGetOwnershipRange( *matrix, &start, &end );

        // for testing
        //std::cerr << "\ngrid size: " << gridSize_ << ", ownership range: [" << start << ", " << end << "]\n\n";

        /* 
         * assemble the matrix elements
         */
        PetscScalar values[ 3 ];
        PetscInt colIndices[ 3 ];
        // This is propably not very efficient...
        for ( PetscInt localI = start; localI < end; ++localI ) {
            if ( localI == 0 ) {
                // the first block
                values[ 0 ] = 1.;
                colIndices[ 0 ] = 0;
                Dune::Petsc::MatSetValues( *matrix, 1, &localI, 1, colIndices, values, INSERT_VALUES );
            } else if ( localI == static_cast< PetscInt >( gridSize_ ) - 1 ) {
                // the last block
                values[ 0 ] = 1.;
                colIndices[ 0 ] = gridSize_ - 1;
                Dune::Petsc::MatSetValues( *matrix, 1, &localI, 1, colIndices, values, INSERT_VALUES );
            } else {
                // a middle block
                values[ 0 ] = -1.;  values[ 1 ] = 2.;  values[ 2 ] = -1.;
                colIndices[ 0 ] = localI - 1;
                colIndices[ 1 ] = localI;
                colIndices[ 2 ] = localI + 1;
                Dune::Petsc::MatSetValues( *matrix, 1, &localI, 3, colIndices, values, INSERT_VALUES );
            }
        }

        // MatAssemblyBegin should be called _after_ all MatSetValues calls...
        Dune::Petsc::MatAssemblyBegin( *matrix, MAT_FINAL_ASSEMBLY );
        Dune::Petsc::MatAssemblyEnd( *matrix, MAT_FINAL_ASSEMBLY );
    }

    template< typename S >
    void assembleRHSSertial ( const Function< S > &f, Vec *vector ) const {
        assert( f.lowerBound() == a_ && f.upperBound() == b_ );

        Dune::Petsc::VecCreate( vector );
        Dune::Petsc::VecSetSizes( *vector, PETSC_DECIDE, static_cast< PetscInt >( gridSize_ ) );
        Dune::Petsc::VecSetFromOptions( *vector );

        // set first and last entry
        Dune::Petsc::VecSetValue( *vector, 0, static_cast< PetscScalar >( 0 ), INSERT_VALUES );
        Dune::Petsc::VecSetValue( *vector, gridSize_ - 1, static_cast< PetscScalar >( 0 ), INSERT_VALUES );

        // set the other entries
        PetscScalar currentX = a_;
        for ( unsigned int i=1; i < gridSize_ - 1; ++i ) {
            currentX += h_;
            Dune::Petsc::VecSetValue( *vector, i, hSquare_ * static_cast< PetscScalar >( f.evaluate( currentX ) ), INSERT_VALUES );
        }

        Dune::Petsc::VecAssemblyBegin( *vector );
        Dune::Petsc::VecAssemblyEnd( *vector );
    }

    template< typename S >
    void assembleRHSParallel ( const Function< S > &f, Vec *vector ) const {
        assert( f.lowerBound() == a_ && f.upperBound() == b_ );

        Dune::Petsc::VecCreate( vector );
        Dune::Petsc::VecSetSizes( *vector, PETSC_DECIDE, static_cast< PetscInt >( gridSize_ ) );
        Dune::Petsc::VecSetFromOptions( *vector );
        
        PetscInt start, end;
        Dune::Petsc::VecGetOwnershipRange( *vector, &start, &end );
        PetscScalar currentX = a_ + start*h_;
        for ( PetscInt localI = start; localI < end; ++localI ) {
            Dune::Petsc::VecSetValue( *vector, localI, hSquare_ * static_cast< PetscScalar >( f.evaluate( currentX ) ), INSERT_VALUES );
            currentX += h_;
        }

        Dune::Petsc::VecAssemblyBegin( *vector );
        Dune::Petsc::VecAssemblyEnd( *vector );
    }

    void initializeSolutionVectorSerial ( Vec *vector ) const {
        Dune::Petsc::VecCreate( vector );
        Dune::Petsc::VecSetSizes( *vector, PETSC_DECIDE, static_cast< PetscInt >( gridSize_ ) );
        Dune::Petsc::VecSetFromOptions( *vector );
        Dune::Petsc::VecSet( *vector, 0. );
    }

    void initializeSolutionVectorParallel ( Vec *vector ) const {
        Dune::Petsc::VecCreate( vector );
        Dune::Petsc::VecSetSizes( *vector, PETSC_DECIDE, static_cast< PetscInt >( gridSize_ ) );
        Dune::Petsc::VecSetFromOptions( *vector );
        Dune::Petsc::VecSet( *vector, 0. );
    }


    /*
     * data fields
     */


    // the number of parts [a,b] is divided into
    unsigned int gridSize_;
    PetscScalar a_;
    PetscScalar b_;
    PetscScalar h_;
    PetscScalar hSquare_;
    bool multipleProcesses_;
};

#endif // #if defined HAVE_PETSC

#endif // LAPLACEASSEMBLER_HH
