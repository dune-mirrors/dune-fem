// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCCOMMON_HH
#define DUNE_FEM_PETSCCOMMON_HH

/*
 * This file should be the first include in every dune-fem-petsc related header
 * It contains some common code in order to deal with PETSc in a consistent, object oriented
 * fashion
 */

#include <string>
#include <iostream>

#include <dune/common/exceptions.hh>

#if HAVE_PETSC

  /* 
   * This turns off PETSc logging of MPI calls. If it is on, PETSc redifines MPI functions as macros, 
   * which breaks some code. E.g. MPI_Allreduce is redefined
   */
#define PETSC_HAVE_BROKEN_RECURSIVE_MACRO

#include <petsc.h>

namespace Dune
{

  namespace Petsc
  {
    // PETSc #defines VecType to be char* - this can cause problems with e.g. ALUGrid 1.50. So we #undef 
    // it here and use a typedef instead. The following trick is really nasty, but conforming to the standard and
    // - most important of all - working :)
    typedef VecType
#undef VecType
    VecType;


#if HAVE_MPI 
#define FEM_PETSC_COMM_DEFAULT  PETSC_COMM_WORLD 
#else
#define FEM_PETSC_COMM_DEFAULT  PETSC_COMM_SELF 
#endif

    /*
     * exceptions
     */
    class Exception : public Dune::Exception {};

    /*
     * types
     */
    typedef ::PetscErrorCode ErrorCode;



    /* The following 2 methods are needed to make the code work with the CHKERRQ
     * macro (which we want to use because it can print useful diagnostic output in
     * case of errors) 
     */
    inline ErrorCode ErrorCheckHelper ( ErrorCode errorCode ) { CHKERRQ( errorCode ); return 0; }

    inline ErrorCode ErrorHandler ( MPI_Comm comm, int line, const char *function, const char *file, const char *dir, 
                                    ErrorCode errorCode, PetscErrorType p, const char *message, void *context )
    {
      std::ostringstream msgout;
      msgout << "PETSc Error in the PETSc function '" << function << "' at " << file << ":" << line << ":";

      const char *text;
      PetscErrorMessage( errorCode, &text, PETSC_NULL );
      if( text )
        msgout << " '" << text << "'. Error message: '" << message << "'";
      else
        msgout << " Unable to retrieve error text";

      DUNE_THROW( Exception, msgout.str() );

      return errorCode;
    }

    inline void ErrorCheck ( ErrorCode errorCode )
    {
      if( errorCode )
      {
        DUNE_THROW( Exception, "Generic PETSC exception" );
      }
    }

    /*
     * This should be called right after the initialization of the MPI manager. It expects the same arguments
     * as PetscInitialize
     */
    inline void initialize( const bool verbose, int &argc, char **&args, const char file[] = 0 , const char help[] = 0, bool ownHandler = true ) 
    {
      ::PetscInitialize( &argc, &args, file, help );

      if( ownHandler )
      {
        // set up our error handler
        if( verbose )
        {
          dverb << "INFORMATION: Setting up an own error handler for PETSc errors. If you want the default PETSc handler,\n"
                << "INFORMATION: set the last argument of Dune::Petsc::initialize(...) to 'false'.\n";
        }
        ::PetscPushErrorHandler( &ErrorHandler, 0 );
      }
    }

#if 0
    /*
     * This should be called right after the initialization of the MPI manager. It expects the same arguments
     * as PetscInitialize
     */
    inline void initialize( int *argc, char ***args, const char file[], const char help[], bool ownHandler = true ) 
    {
      ::PetscInitialize( argc, args, file, help );

      if( ownHandler )
      {
        // set up our error handler
        if( Dune::Fem::MPIManager::rank() == 0 )
        {
          std::cerr << "INFORMATION: Setting up an own error handler for PETSc errors. If you want the default PETSc handler,\n"
                    << "INFORMATION: set the last argument of Dune::Petsc::initialize(...) to 'false'.\n";
        }
        ::PetscPushErrorHandler( &ErrorHandler, 0 );
      }
    }
#endif

    /*
     * This should be called just before the termination of the program.
     */
    inline void finalize ()
    {
      // TODO: test here if we are using our own handler
      ::PetscPopErrorHandler();

      ::PetscFinalize();
    }

    /*
     * ==================================================
     * These are simple mappings to PETSc's C-routines. (Maybe some of them are not needed...).
     *
     * The PETSC_VERSION_... customizations are not very well tested yet
     */
    inline void KSPCreate ( KSP *inksp ) { ErrorCheck( ::KSPCreate( FEM_PETSC_COMM_DEFAULT, inksp ) ); }
    inline void KSPDestroy ( KSP *ksp ) { 
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::KSPDestroy( *ksp ) ); 
#else
        ErrorCheck( ::KSPDestroy( ksp ) ); 
#endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    inline void KSPGetPC ( KSP ksp, PC *pc ) { ErrorCheck( ::KSPGetPC( ksp, pc ) ); }
    inline void KSPSetFromOptions ( KSP ksp ) { ErrorCheck( ::KSPSetFromOptions( ksp ) ); }
    inline void KSPSetType ( KSP ksp, const KSPType type ) { ErrorCheck( ::KSPSetType( ksp, type ) ); }
    inline void KSPGMRESSetRestart ( KSP ksp, PetscInt restart ) { ErrorCheck( ::KSPGMRESSetRestart( ksp, restart ) ); }
    inline void KSPView ( KSP ksp, PetscViewer viewer = PETSC_VIEWER_STDOUT_(FEM_PETSC_COMM_DEFAULT)  ) 
    { 
      ErrorCheck( ::KSPView( ksp, viewer ) ); 
    }
    inline void KSPMonitorSet (KSP ksp, PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
                               void *mctx,PetscErrorCode (*monitordestroy)(void*) 
#else
                               void *mctx,PetscErrorCode (*monitordestroy)(void**)  
#endif
                       )
    { 
        ErrorCheck( ::KSPMonitorSet( ksp, monitor, mctx, monitordestroy ) ); 
    }
    inline void KSPGetIterationNumber( KSP ksp, PetscInt* its ) 
    { ErrorCheck( ::KSPGetIterationNumber( ksp, its ) ); }

    inline void KSPSetOperators (KSP ksp, Mat Amat, Mat Pmat, MatStructure flag ) { ErrorCheck( ::KSPSetOperators( ksp, Amat, Pmat, flag ) ); }
    inline void KSPSetTolerances ( KSP ksp, PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits ) 
      { ErrorCheck( ::KSPSetTolerances( ksp, rtol, abstol, dtol, maxits ) ); }
    inline void KSPSolve ( KSP ksp, Vec b, Vec x ) { ErrorCheck( ::KSPSolve( ksp, b, x ) ); }
    inline void KSPSetPC ( KSP ksp, PC pc ) { ErrorCheck( ::KSPSetPC( ksp, pc ) ); }

    // preconditioning 
    inline void PCCreate  ( PC* pc) { ErrorCheck( ::PCCreate( FEM_PETSC_COMM_DEFAULT, pc ) ); }
    inline void PCDestroy ( PC* pc) { 
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
      ErrorCheck( ::PCDestroy( *pc ) ); 
#else  
      ErrorCheck( ::PCDestroy( pc ) ); 
#endif 
    }
    inline void PCSetType ( PC pc, const PCType type ) { ErrorCheck( ::PCSetType(  pc, type ) ); }
    inline void PCFactorSetLevels( PC pc, PetscInt level ) { ErrorCheck( ::PCFactorSetLevels(  pc, level ) ); }

    // matrix routines 
    inline void MatAssemblyBegin ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyBegin( mat, type ) ); }
    inline void MatAssemblyEnd ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyEnd( mat, type ) ); }
    inline void MatCreate ( Mat *A ) { ErrorCheck( ::MatCreate( FEM_PETSC_COMM_DEFAULT, A) ); }
    inline void MatCreateBlockMat ( Mat *A, PetscInt m, PetscInt n, PetscInt bs, PetscInt nz, PetscInt* nnz ) 
    { 
      ErrorCheck( ::MatCreateBlockMat( FEM_PETSC_COMM_DEFAULT, n, m, bs, nz, nnz, A) ); 
    }
    inline void MatDestroy ( Mat *A ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::MatDestroy( *A ) ); 
      #else
        ErrorCheck( ::MatDestroy( A ) ); 
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    inline void MatSetUp( Mat mat ) 
    { 
      ErrorCheck( ::MatSetUp(mat)); 
    }
    inline void MatSetUp( Mat mat, PetscInt bs, int nz ) 
    { 
      if (bs == 1)
      {
        ErrorCheck( ::MatSeqAIJSetPreallocation(mat,nz,PETSC_NULL) );
        ErrorCheck( ::MatMPIAIJSetPreallocation(mat,nz,PETSC_NULL,nz/2,PETSC_NULL) );
      }
      else
      {
        ErrorCheck( ::MatSeqBAIJSetPreallocation(mat,bs,nz,PETSC_NULL) );
        ErrorCheck( ::MatMPIBAIJSetPreallocation(mat,bs,nz,PETSC_NULL,nz/2,PETSC_NULL) );
      }
      // the following seems not to work for block matrix
      ErrorCheck( ::MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) );
      // the following only works for block matrix
      // but should be used with MAT_NEW_NONZERO_LOCATIONS...
      // ErrorCheck( ::MatSetOption(mat, MAT_USE_HASH_TABLE,PETSC_FALSE) );
      ErrorCheck( ::MatSetUp(mat)); 
    }
    inline void MatSetUp( Mat mat, PetscInt bs, const int *d_nnz, const int *o_nnz )
    {
      if (bs == 1)
      {
        ErrorCheck( ::MatSeqAIJSetPreallocation(mat,0,d_nnz ) );
        ErrorCheck( ::MatMPIAIJSetPreallocation(mat,0,d_nnz,5,o_nnz) );
      }
      else
      {
        ErrorCheck( ::MatSeqBAIJSetPreallocation(mat,bs,0,d_nnz ) );
        ErrorCheck( ::MatMPIBAIJSetPreallocation(mat,bs,0,d_nnz,5,PETSC_NULL) );
      }
      // see previous comments
      ErrorCheck( ::MatSetOption(mat, MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_FALSE) );
      ErrorCheck( ::MatSetUp(mat));
    }

    inline void MatGetOwnershipRange ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetOwnershipRange( mat, m, n ) ); }
    inline void MatGetSize ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetSize( mat, m, n ) ); }
    inline void MatMult  ( Mat mat, Vec x, Vec y ) { ErrorCheck( ::MatMult( mat, x, y ) ); }
    inline void MatSetBlockSize ( Mat A, PetscInt bs ) { ErrorCheck( ::MatSetBlockSize( A, bs ) ); }
    inline void MatSetSizes ( Mat A, PetscInt m, PetscInt n, PetscInt M, PetscInt N ) { ErrorCheck( ::MatSetSizes( A, m, n, M, N ) ); }
    inline void MatSetFromOptions ( Mat B ) { ErrorCheck( ::MatSetFromOptions( B ) ); }
    inline void MatSetType ( Mat mat, const MatType matype ) { ErrorCheck( ::MatSetType( mat, matype ) ); }
    inline void MatSetValue ( Mat v, PetscInt i, PetscInt j, PetscScalar va, InsertMode mode ) { ErrorCheck( ::MatSetValue( v, i, j, va, mode ) ); }
    inline void MatSetValues ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv )
      { ErrorCheck( ::MatSetValues( mat, m, idxm, n, idxn, v, addv ) ); }
    inline void MatGetValues ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[] )
      { ErrorCheck( ::MatGetValues( mat, m, idxm, n, idxn, v ) ); }
    inline void MatView ( Mat mat, PetscViewer viewer ) { ErrorCheck( ::MatView( mat, viewer ) ); }
    inline void MatZeroEntries ( Mat mat ) { ErrorCheck( ::MatZeroEntries( mat ) ); }
    inline void PetscBarrier ( PetscObject obj ) { ErrorCheck( ::PetscBarrier( obj ) ); }
    inline void PetscFinalize () { ErrorCheck( ::PetscFinalize() ); }
    inline void PetscInitialize( int *argc, char ***args, const char file[], const char help[] ) { ErrorCheck( ::PetscInitialize( argc, args, file, help ) ); }
    inline void PetscViewerASCIIOpen ( MPI_Comm comm, const char name[], PetscViewer *lab ) { ErrorCheck( ::PetscViewerASCIIOpen( comm, name, lab ) ); }
    inline void PetscViewerDestroy ( PetscViewer *viewer ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::PetscViewerDestroy( *viewer ) );  
      #else
        ErrorCheck( ::PetscViewerDestroy( viewer ) );  
      #endif
    }
    inline void PetscViewerSetFormat ( PetscViewer viewer, PetscViewerFormat format ) { ErrorCheck( ::PetscViewerSetFormat( viewer, format ) ); }
    inline void VecAssemblyBegin ( Vec vec ) { ErrorCheck( ::VecAssemblyBegin( vec ) ); }
    inline void VecAssemblyEnd ( Vec vec ) { ErrorCheck( ::VecAssemblyEnd( vec ) ); }
    inline void VecAXPY ( Vec y, PetscScalar alpha, Vec x) { ErrorCheck( ::VecAXPY( y, alpha, x ) ); }
    inline void VecCopy ( Vec x, Vec y ) { ErrorCheck( ::VecCopy( x, y ) ); }
    inline void VecCreate ( Vec *vec ) { ErrorCheck( ::VecCreate( FEM_PETSC_COMM_DEFAULT, vec ) ); }
    inline void VecCreateGhost ( PetscInt n, PetscInt N, PetscInt nghost, const PetscInt ghosts[], Vec *vv )
      { ErrorCheck( ::VecCreateGhost( FEM_PETSC_COMM_DEFAULT, n, N, nghost, ghosts, vv ) ); }
    inline void VecDestroy ( Vec *v ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::VecDestroy( *v ) ); 
      #else
        ErrorCheck( ::VecDestroy( v ) ); 
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    } 
    inline void VecDot ( Vec x, Vec y, PetscScalar *val ) { ErrorCheck( ::VecDot( x, y, val ) ); }
    inline void VecDuplicate ( Vec v, Vec *newv ) { ErrorCheck( ::VecDuplicate( v, newv ) ); }
    inline void VecGetLocalSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetLocalSize( x, size ) ); }
    inline void VecGetOwnershipRange ( Vec x, PetscInt *low, PetscInt *high ) { ErrorCheck( ::VecGetOwnershipRange( x, low, high ) ); }
    inline void VecGetSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetSize( x, size ) ); }
    inline void VecGetValues ( Vec x, PetscInt ni, const PetscInt ix[], PetscScalar y[] ) { ErrorCheck( ::VecGetValues( x, ni, ix, y ) ); }
    inline void VecGhostGetLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostGetLocalForm( g, l ) ); }
    inline void VecGhostRestoreLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostRestoreLocalForm( g, l ) ); }
    inline void VecGhostUpdateBegin ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateBegin( g, insertmode, scattermode ) ); }
    inline void VecGhostUpdateEnd ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateEnd( g, insertmode, scattermode ) ); }
    inline void VecNorm ( Vec x, NormType type, PetscReal *val ) { ErrorCheck( ::VecNorm( x, type, val ) ); }
    inline void VecScale ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecScale( x, alpha ) ); }
    inline void VecSet ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecSet( x, alpha ) ); }
    inline void VecSetFromOptions ( Vec vec ) { ErrorCheck( ::VecSetFromOptions( vec ) ); }
    inline void VecSetType ( Vec vec, const VecType method ) { ErrorCheck( ::VecSetType( vec, method ) ); }
    inline void VecSetSizes ( Vec v, PetscInt n, PetscInt N ) { ErrorCheck( ::VecSetSizes( v, n, N ) ); }
    inline void VecSetValue ( Vec v, int row, PetscScalar value, InsertMode mode ) { ErrorCheck( ::VecSetValue( v, row, value, mode ) ); }
    inline void VecView ( Vec vec, PetscViewer viewer ) { ErrorCheck( ::VecView( vec, viewer ) ); }
      
  } // namespace Petsc

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCCOMMON_HH
