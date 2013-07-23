#ifndef L2ERROR_HH
#define L2ERROR_HH

#include <assert.h>
#include <iostream>
#include <sstream>

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

#include "function.hh"
#include "printer.hh"

// calculate the l2 norm (or some norm equivalent to this, I dunno) of numericalSolution-exactSolution
template< typename S >
class L2Error {
public:
    typedef S ScalarType;

    template< typename T >
    static ScalarType calculateError ( Vec numericalSolution, const Function< T > &exactSolution, 
                                       unsigned int gridSize, int eocStep = 0 )
    {
        #ifndef NDEBUG
        PetscInt vecLength;
        Dune::Petsc::VecGetSize( numericalSolution, &vecLength );
        assert( static_cast< unsigned int >( vecLength ) == gridSize );
        #endif // NDEBUG

        Vec newVec;
        Dune::Petsc::VecDuplicate( numericalSolution, &newVec );
        Dune::Petsc::VecCopy( numericalSolution, newVec );

        ScalarType h = ( exactSolution.upperBound() - exactSolution.lowerBound() ) / static_cast< ScalarType >( gridSize - 1 );
        ScalarType a = exactSolution.lowerBound();
        
        PetscInt start, end;
        Dune::Petsc::VecGetOwnershipRange( newVec, &start, &end );
        //std::cerr << "owner-range in calculateError: start: " << start << ", end: " << end << std::endl;
    
        ScalarType currentX = a + start*h;
        for ( PetscInt localI = start; localI < static_cast< PetscInt >( end ); ++localI ) {
            Dune::Petsc::VecSetValue( newVec, localI, -exactSolution.evaluate( currentX ), ADD_VALUES );
            currentX += h;
        }
           
        Dune::Petsc::VecAssemblyBegin( newVec );
        Dune::Petsc::VecAssemblyEnd( newVec );

        std::stringstream gnuplotFilename;
        gnuplotFilename << "error.out.gnuplot." << eocStep;
        std::ofstream outFile( gnuplotFilename.str() );
        Printer::printVec( gridSize, exactSolution.lowerBound(), exactSolution.upperBound(), newVec, outFile );
        

        PetscReal ret = 0;
        Dune::Petsc::VecNorm( newVec, NORM_2, &ret );
        ret *= std::sqrt( h );

        return static_cast< ScalarType >( ret );
    }

private:
    L2Error ();
};

#endif // #if defined HAVE_PETSC

#endif // L2ERROR_HH
