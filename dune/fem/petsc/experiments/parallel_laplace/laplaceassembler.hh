#ifndef LAPLACEASSEMBLER_HH
#define LAPLACEASSEMBLER_HH

#include <assert.h>
#include <iostream>
#include <sstream>
#include <algorithm>

#if defined HAVE_PETSC

#include <dune/fem/petsc/common/petsccommon.hh>

#include "petscmat.h"
#include "petscvec.h"

#include "function.hh"
#include "printer.hh"


/*
 * Class for PETSc matrix and vector assembly
 */
//template< bool multipleProcesses >
class LaplaceAssembler {
public:

    LaplaceAssembler ( unsigned int gridSize, PetscScalar a, PetscScalar b )
    : gridSize_( gridSize ),
      a_( a ),
      b_( b ),
      h_( static_cast< PetscScalar >( b_ - a_ )/ static_cast< PetscScalar >( gridSize_ - 1 ) ),
      hSquare_( h_ * h_ )
    {
        assert( b_ > a_ );
    }

    // Creates and assembles the matrix 'matrix'. Calls, among other functions, MatCreate(...) on
    // 'matrix', so 'matrix' should be supplied entirely uninitialized...
    void assembleMatrix ( Mat *matrix ) const {

        int myRank;
        int commSize;
        MPI_Comm_rank( PETSC_COMM_WORLD, &myRank );
        MPI_Comm_size( PETSC_COMM_WORLD, &commSize );

        int myNumDofs = getVecLocalSize( commSize, myRank );

        int startIndex = 0;
        MPI_Scan( &myNumDofs, &startIndex, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD );
        startIndex -= myNumDofs;

        Dune::Petsc::MatCreate( matrix );
        Dune::Petsc::MatSetSizes( *matrix, myNumDofs, myNumDofs, PETSC_DECIDE, PETSC_DECIDE );
        Dune::Petsc::MatSetType( *matrix, MATMPIAIJ );

        PetscInt glRows;
        MPI_Allreduce( &myNumDofs, &glRows, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD );

        int limit = startIndex + myNumDofs;

        PetscScalar firstRow[ 2 ];
        PetscScalar secondRow[ 2 ];
        PetscScalar one = 1.;
        PetscInt colIndices[ 2 ];
        firstRow[ 0 ] = 1.; firstRow[ 1 ] = -1.;
        secondRow[ 0 ] = -1.; secondRow[ 1 ] = 1.;
        PetscInt i;
        PetscInt iPlusOne;

        if( startIndex == 0 ) {
            // so we are on the first process
            startIndex++;
            i = 0;
            iPlusOne = i+1;
            colIndices[ 0 ] = i; colIndices[ 1 ] = i+1;
            Dune::Petsc::MatSetValues( *matrix, 1, &i, 1, &i, &one, ADD_VALUES );
            Dune::Petsc::MatSetValues( *matrix, 1, &iPlusOne, 2, colIndices, secondRow, ADD_VALUES );
        }

        if( limit == glRows ) {
            // so we are on the last process
            limit -= 2;
            i = limit;
            iPlusOne = i+1;
            colIndices[ 0 ] = i; colIndices[ 1 ] = i+1;
            Dune::Petsc::MatSetValues( *matrix, 1, &i, 2, colIndices, firstRow, ADD_VALUES );
            Dune::Petsc::MatSetValues( *matrix, 1, &iPlusOne, 1, &iPlusOne, &one, ADD_VALUES );
        }

        for( i = startIndex; i < limit; ++i ) {
            iPlusOne = i+1;
            colIndices[ 0 ] = i; colIndices[ 1 ] = i+1;
            Dune::Petsc::MatSetValues( *matrix, 1, &i, 2, colIndices, firstRow, ADD_VALUES );
            Dune::Petsc::MatSetValues( *matrix, 1, &iPlusOne, 2, colIndices, secondRow, ADD_VALUES );
        }
        
        // MatAssemblyBegin should be called _after_ all MatSetValues calls...
        Dune::Petsc::MatAssemblyBegin( *matrix, MAT_FINAL_ASSEMBLY );
        Dune::Petsc::MatAssemblyEnd( *matrix, MAT_FINAL_ASSEMBLY );

        // Don't know why, but this sais: "Matrix Object: 1 MPI processes", although we use more than 1 process.
        // After stepping through PETSc code I think this is not an error by me. I believe that the Mat is converted to 
        // a sequential AIJ matrix internally... (?)
        //Dune::Petsc::MatView( *matrix, PETSC_VIEWER_STDOUT_WORLD );
    
    }

    // assemble the right hand side vector
    template< typename S >
    void assembleRHS ( const Function< S > &f, Vec *vector, int eocStep = 0 ) const {
        /*
         * We want to assemble the RHS using a ghosted Vec here.
         */
        assert( f.lowerBound() == a_ && f.upperBound() == b_ );

        int myRank;
        int commSize;
        MPI_Comm_rank( PETSC_COMM_WORLD, &myRank );
        MPI_Comm_size( PETSC_COMM_WORLD, &commSize );

        int myNumDofs = getVecLocalSize( commSize, myRank );

        int startIndex;
        MPI_Scan( &myNumDofs, &startIndex, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD );
        startIndex -= myNumDofs;

        // Setup the ghost array
        PetscInt nGhosts;
        PetscInt ghost;
        if( myRank == 0 ) {
            nGhosts = 0;
        } else {
            nGhosts = 1;
            ghost = startIndex - 1;
        }

        Dune::Petsc::VecCreateGhost( myNumDofs, PETSC_DECIDE, nGhosts, &ghost, vector );

        Vec ghostedVec;
        Dune::Petsc::VecGhostGetLocalForm( *vector, &ghostedVec );
        Dune::Petsc::VecSet( ghostedVec, 0. );

        // add to the ghost dof if necessary
        if( myRank > 0 ) {
            PetscScalar val = .5*hSquare_*static_cast< PetscScalar >( f.evaluate(  a_ + ( startIndex - 1 + .5 )*h_ ) );
            Dune::Petsc::VecSetValue( ghostedVec, myNumDofs, val, ADD_VALUES );
            Dune::Petsc::VecSetValue( ghostedVec, 0, val, ADD_VALUES );
        }

        for( PetscInt i = 0; i < myNumDofs - 1; ++i ) {
            PetscScalar val = .5*hSquare_*static_cast< PetscScalar >( f.evaluate( a_ + (startIndex + i + .5 )*h_ ) );
            Dune::Petsc::VecSetValue( ghostedVec, i, val, ADD_VALUES );
            Dune::Petsc::VecSetValue( ghostedVec, i + 1, val, ADD_VALUES );
        }
        Dune::Petsc::VecGhostUpdateBegin( *vector, ADD_VALUES, SCATTER_REVERSE );
        Dune::Petsc::VecGhostUpdateEnd( *vector, ADD_VALUES, SCATTER_REVERSE );
        Dune::Petsc::VecGhostUpdateBegin( *vector, INSERT_VALUES, SCATTER_FORWARD );
        Dune::Petsc::VecGhostUpdateEnd( *vector, INSERT_VALUES, SCATTER_FORWARD );
        Dune::Petsc::VecGhostRestoreLocalForm( *vector, &ghostedVec );

        std::stringstream gnuplotFilename;
        gnuplotFilename << "rhs.out.gnuplot." << eocStep;
        std::ofstream outFile( gnuplotFilename.str() );
        Printer::printVec( gridSize_, f.lowerBound(), f.upperBound(), *vector, outFile );
    }

    void finalizeRHS( Vec *vector ) const {
        Dune::Petsc::VecDestroy( vector );
    }

    void initializeSolutionVector ( Vec *vector ) const {
        int myNumDofs = getVecLocalSize();

        Dune::Petsc::VecCreate( vector );
        Dune::Petsc::VecSetType( *vector, VECMPI );
        Dune::Petsc::VecSetSizes( *vector, myNumDofs, PETSC_DECIDE );
        Dune::Petsc::VecSet( *vector, 0. );
    }

private:
    /*
     * private methods
     */

    // get the number of dofs that a PETSc Vec should own on this process
    int getVecLocalSize( int commSize, int commRank ) const {
        // Determine the local sizes for the vec
        int sz = gridSize_ / commSize;
        int rest = gridSize_ % commSize;
        int myNumDofs = sz;
        // rank 0 has the most dofs...
        if( commRank == 0 ) {
            myNumDofs += rest;
        }

        // assertion for debugging
        int totalSize;
        MPI_Allreduce( &myNumDofs, &totalSize, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD );
        assert( totalSize == int( gridSize_ ) );

        return myNumDofs;
    }

    // same as the version with 2 arguments
    int getVecLocalSize () const {
        int commRank;
        int commSize;
        MPI_Comm_rank( PETSC_COMM_WORLD, &commRank );
        MPI_Comm_size( PETSC_COMM_WORLD, &commSize );
        return getVecLocalSize( commSize, commRank );
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
};

#endif // #if defined HAVE_PETSC

#endif // LAPLACEASSEMBLER_HH
