#ifndef PRINTER_HH
#define PRINTER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdio>

#if defined HAVE_PETSC

#include <dune/fem/misc/petsc/petsccommon.hh>

#include "petscvec.h"
#include "petscmat.h"
#include "petscviewer.h"


#include "function.hh"


class Printer {
public:
    
    // prints both solutions to std::cout in a format that gnuplot understands
    template< typename Stream >
    static void printVec ( unsigned int gridSize, PetscScalar a, PetscScalar b, Vec &vec, Stream &stream )
    {
        const char* bufferFileName = "buffer.out.tmp.5f5285a5192800406c203d678927718c";
        PetscViewer viewer;
        MPI_Comm comm = PETSC_COMM_WORLD;
        Dune::Petsc::PetscViewerASCIIOpen( comm, bufferFileName, &viewer );
        Dune::Petsc::PetscViewerSetFormat ( viewer, PETSC_VIEWER_ASCII_INFO );
        Dune::Petsc::VecView( vec, viewer );
        Dune::Petsc::PetscViewerDestroy( &viewer );

        int myRank;
        MPI_Comm_rank( PETSC_COMM_WORLD, &myRank );
        if( myRank != 0 ) {
            return;
        }

        /* 
         * Read the whole file
         */
        std::vector< std::string > lines;
        std::ifstream vecFile( bufferFileName );
        const int lineSz = 1024;
        char line[ lineSz ];

        // clear the first 3 lines, as they don't contain values
        vecFile.getline( line, lineSz ); vecFile.getline( line, lineSz ); vecFile.getline( line, lineSz );

        do {
            vecFile.getline( line, lineSz );
            // skip the lines containing "Process...."
            if( line[ 0 ] == 'P' ) {
                vecFile.getline( line, lineSz );
            }
            if( vecFile.gcount() > 0 ) {
                lines.push_back( line );
            }
        } while( !vecFile.eof() );

        /*
         * delete the tmp file and write everything to the stream then
         */
        std::remove( bufferFileName );

        assert( lines.size() == gridSize );
        assert( b > a );
        PetscScalar currentX = a;
        PetscScalar h = ( b - a )/gridSize;

        for( int i = 0; i < int( gridSize ); ++i ) {
            stream << currentX << "\t\t" << lines[ i ] << std::endl;
            currentX += h;
        }

    }

private:
    // prohibited methods
    Printer ();

};

#endif // #if defined HAVE_PETSC

#endif // PRINTER_HH
