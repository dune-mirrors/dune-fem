#ifndef DUNE_FEM_PETSCCOMMON_HH
#define DUNE_FEM_PETSCCOMMON_HH

/*
 * This file should be the first include in every dune-fem-petsc related header
 * It contains some common code in order to deal with PETSc in a consistent, object oriented
 * fashion
 */

#include <string>
#include <iostream>
#include <iomanip>

#include <dune/common/parallel/communication.hh>
#include <dune/common/stdstreams.hh>
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
    inline static ErrorCode ErrorCheckHelper ( ErrorCode errorCode ) { CHKERRQ( errorCode ); return 0; }

    inline ErrorCode ErrorHandler ( MPI_Comm comm, int line, const char *function, const char *file,
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
                                    const char *dir,
#endif
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

    inline static void ErrorCheck ( ErrorCode errorCode )
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
    inline bool initialize( const bool verbose, int &argc, char **&args, const char file[] = 0 , const char help[] = 0, bool ownHandler = true )
    {
      bool wasInitializedHere = false ;
      PetscBool petscInitialized = PETSC_FALSE;

      // check whether PETSc had been initialized elsewhere
      ::PetscInitialized( &petscInitialized );

      if( ! petscInitialized )
      {
        ::PetscInitialize( &argc, &args, file, help );
        wasInitializedHere = true;
      }

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
      return wasInitializedHere;
    }

    /*
     * This should be called just before the termination of the program.
     */
    inline void finalize ()
    {
      // TODO: test here if we are using our own handler
      ::PetscPopErrorHandler();
      PetscBool finalized = PETSC_FALSE;
      ErrorCheck( ::PetscFinalized( &finalized ) );

      if( ! finalized )
      {
        ::PetscFinalize();
      }
    }

    template <class Comm>
    inline auto castToPetscComm( const Comm& comm )
    {
      // this is needed because Dune::No_Comm or
      // Dune::Communication< No_Comm > does not cast into MPI_Comm
      if constexpr ( std::is_same< Dune::No_Comm, Comm > :: value ||
                     std::is_same< Dune::Communication< No_Comm >, Comm >::value )
      {
        return PETSC_COMM_SELF;
      }
      else
      {
        return comm;
      }
    }

    /*
     * ==================================================
     * These are simple mappings to PETSc's C-routines. (Maybe some of them are not needed...).
     *
     * The PETSC_VERSION_... customizations are not very well tested yet
     */
    template <class Comm>
    inline static void KSPCreate ( const Comm& comm, KSP *inksp )
    {
      ErrorCheck( ::KSPCreate( castToPetscComm( comm ), inksp ) );
    }

    inline static void KSPDestroy ( KSP *ksp ) {
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::KSPDestroy( *ksp ) );
#else
        ErrorCheck( ::KSPDestroy( ksp ) );
#endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    inline static void KSPGetPC ( KSP ksp, PC *pc ) { ErrorCheck( ::KSPGetPC( ksp, pc ) ); }
    inline static void KSPSetFromOptions ( KSP ksp ) { ErrorCheck( ::KSPSetFromOptions( ksp ) ); }
    inline static void KSPSetUp ( KSP ksp ) { ErrorCheck( ::KSPSetUp( ksp ) ); }
    inline static void KSPSetType ( KSP ksp, const KSPType type ) { ErrorCheck( ::KSPSetType( ksp, type ) ); }
    inline static void KSPGMRESSetRestart ( KSP ksp, PetscInt restart ) { ErrorCheck( ::KSPGMRESSetRestart( ksp, restart ) ); }

    template <class Comm>
    inline static void KSPView ( const Comm& comm,
                                 KSP ksp,
                                 PetscViewer viewer )
    {
      ErrorCheck( ::KSPView( ksp, viewer ) );
    }

    template <class Comm>
    inline static void KSPView ( const Comm& comm,
                                 KSP ksp )
    {
      KSPView( comm, ksp, PETSC_VIEWER_STDOUT_( castToPetscComm(comm) ) );
    }

    inline static void KSPMonitorSet (KSP ksp, PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
                               void *mctx,PetscErrorCode (*monitordestroy)(void*)
#else
                               void *mctx,PetscErrorCode (*monitordestroy)(void**)
#endif
                       )
    {
        ErrorCheck( ::KSPMonitorSet( ksp, monitor, mctx, monitordestroy ) );
    }
    inline static void KSPGetIterationNumber( KSP ksp, PetscInt* its )
    { ErrorCheck( ::KSPGetIterationNumber( ksp, its ) ); }
    inline static void KSPGetConvergedReason(KSP ksp, KSPConvergedReason *reason)
    { ErrorCheck( ::KSPGetConvergedReason( ksp, reason ) ); }


#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
    inline static void KSPSetOperators (KSP ksp, Mat Amat, Mat Pmat, MatStructure flag ) { ErrorCheck( ::KSPSetOperators( ksp, Amat, Pmat, flag ) ); }
#else
    inline static void KSPSetOperators (KSP ksp, Mat Amat, Mat Pmat ) { ErrorCheck( ::KSPSetOperators( ksp, Amat, Pmat ) ); }
#endif
    inline static void KSPSetTolerances ( KSP ksp, PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits )
      { ErrorCheck( ::KSPSetTolerances( ksp, rtol, abstol, dtol, maxits ) ); }
    inline static void KSPSetInitialGuessNonzero( KSP ksp, PetscBool flg ) { ErrorCheck( ::KSPSetInitialGuessNonzero( ksp, flg ) ); };
    inline static void KSPSolve ( KSP ksp, Vec b, Vec x ) { ErrorCheck( ::KSPSolve( ksp, b, x ) ); }
    inline static void KSPSetPC ( KSP ksp, PC pc ) { ErrorCheck( ::KSPSetPC( ksp, pc ) ); }

    // preconditioning
    template <class Comm>
    inline static void PCCreate  ( const Comm& comm, PC* pc)
    {
      ErrorCheck( ::PCCreate( castToPetscComm( comm ), pc ) );
    }

    inline static void PCDestroy ( PC* pc) {
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
      ErrorCheck( ::PCDestroy( *pc ) );
#else
      ErrorCheck( ::PCDestroy( pc ) );
#endif
    }
    inline static void PCSetType ( PC pc, const PCType type ) { ErrorCheck( ::PCSetType(  pc, type ) ); }
    inline static void PCSetFromOptions ( PC pc ) { ErrorCheck( ::PCSetFromOptions(  pc ) ); }
    inline static void PCSetUp ( PC pc ) { ErrorCheck( ::PCSetUp(  pc ) ); }
    inline static void PCFactorSetLevels( PC pc, PetscInt level ) { ErrorCheck( ::PCFactorSetLevels(  pc, level ) ); }
    inline static void PCSORSetOmega( PC pc, PetscReal omega ) { ErrorCheck( ::PCSORSetOmega(  pc, omega ) ); }
    inline static void PCSORSetSymmetric( PC pc, MatSORType flag ) { ErrorCheck( ::PCSORSetSymmetric(  pc, flag ) ); }
#if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 9
    inline static void PCFactorSetMatSolverPackage( PC pc, const MatSolverPackage type )
    {
      ErrorCheck( ::PCFactorSetMatSolverPackage(  pc, type ) );
    }
#else
    inline static void PCFactorSetMatSolverPackage( PC pc, const MatSolverType type )
    {
      ErrorCheck( ::PCFactorSetMatSolverType( pc, type ) );
    }
#endif
    inline static void PCHYPRESetType( PC pc, const char* type )
    {
      ErrorCheck( ::PCHYPRESetType( pc, type ) );
    }

    // matrix routines
    inline static void MatAssemblyBegin ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyBegin( mat, type ) ); }
    inline static void MatAssemblyEnd ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyEnd( mat, type ) ); }
    inline static void MatAssembled( Mat mat, PetscBool* assembled ) { ErrorCheck( ::MatAssembled( mat, assembled ) ); }

    template <class Comm>
    inline static void MatCreate ( const Comm& comm, Mat *A )
    {
      ErrorCheck( ::MatCreate( castToPetscComm( comm ), A) );
    }

    template <class Comm>
    inline static void MatCreateBlockMat ( const Comm& comm,
                                           Mat *A,
                                           PetscInt m, PetscInt n, PetscInt bs, PetscInt nz, PetscInt* nnz )
    {
      ErrorCheck( ::MatCreateBlockMat( castToPetscComm( comm ), n, m, bs, nz, nnz, A) );
    }
    inline static void MatDestroy ( Mat *A ) {
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::MatDestroy( *A ) );
      #else
        ErrorCheck( ::MatDestroy( A ) );
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    inline static void MatSetUp( Mat mat )
    {
      ErrorCheck( ::MatSetUp(mat));
    }
    inline static void MatSetUp( Mat mat, PetscInt bs, PetscInt nz )
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
    inline static void MatSetUp( Mat mat, PetscInt bs, const PetscInt *d_nnz, const PetscInt *o_nnz )
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

    inline static void MatGetOwnershipRange ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetOwnershipRange( mat, m, n ) ); }
    inline static void MatGetSize ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetSize( mat, m, n ) ); }
    inline static void MatMult  ( Mat mat, Vec x, Vec y ) { ErrorCheck( ::MatMult( mat, x, y ) ); }
    inline static void MatSetBlockSize ( Mat A, PetscInt bs ) { ErrorCheck( ::MatSetBlockSize( A, bs ) ); }
    inline static void MatSetSizes ( Mat A, PetscInt m, PetscInt n, PetscInt M, PetscInt N ) { ErrorCheck( ::MatSetSizes( A, m, n, M, N ) ); }
    inline static void MatSetFromOptions ( Mat B ) { ErrorCheck( ::MatSetFromOptions( B ) ); }
    inline static void MatSetType ( Mat mat, const MatType matype ) { ErrorCheck( ::MatSetType( mat, matype ) ); }

    inline static void MatSetValue ( Mat v, PetscInt i, PetscInt j, PetscScalar va, InsertMode mode )
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        ErrorCheck( ::MatSetValue( v, i, j, va, mode ) );
      }
    }

    inline static void MatSetValues ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv )
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        ErrorCheck( ::MatSetValues( mat, m, idxm, n, idxn, v, addv ) );
      }
    }

    inline static void MatSetValuesBlocked ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv )
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        ErrorCheck( ::MatSetValuesBlocked( mat, m, idxm, n, idxn, v, addv ) );
      }
    }

    inline static void MatGetValues ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], PetscScalar v[] )
      { ErrorCheck( ::MatGetValues( mat, m, idxm, n, idxn, v ) ); }

    inline static void MatZeroRows ( Mat mat, PetscInt m, const PetscInt idxm[], const PetscScalar v )
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        ErrorCheck( ::MatZeroRows( mat, m, idxm, v, 0, 0 ) );
      }
    }

    inline static void MatView ( Mat mat, PetscViewer viewer ) { ErrorCheck( ::MatView( mat, viewer ) ); }
    inline static void MatZeroEntries ( Mat mat )
    {
#ifdef _OPENMP
#pragma omp critical
#endif
      {
        ErrorCheck( ::MatZeroEntries( mat ) );
      }
    }

    inline static void PetscBarrier ( PetscObject obj ) { ErrorCheck( ::PetscBarrier( obj ) ); }
    inline static void PetscFinalize () { ErrorCheck( ::PetscFinalize() ); }
    inline static void PetscInitialize( int *argc, char ***args, const char file[], const char help[] ) { ErrorCheck( ::PetscInitialize( argc, args, file, help ) ); }
    inline static void PetscViewerASCIIOpen ( MPI_Comm comm, const char name[], PetscViewer *lab ) { ErrorCheck( ::PetscViewerASCIIOpen( comm, name, lab ) ); }
    inline static void PetscViewerDestroy ( PetscViewer *viewer ) {
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::PetscViewerDestroy( *viewer ) );
      #else
        ErrorCheck( ::PetscViewerDestroy( viewer ) );
      #endif
    }
    inline static void PetscViewerSetFormat ( PetscViewer viewer, PetscViewerFormat format )
    {
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 5
        ErrorCheck( ::PetscViewerSetFormat( viewer, format ) );
      #else
        ErrorCheck( ::PetscViewerPushFormat( viewer, format ) );
      #endif
    }
    inline static void VecAssemblyBegin ( Vec vec ) { ErrorCheck( ::VecAssemblyBegin( vec ) ); }
    inline static void VecAssemblyEnd ( Vec vec ) { ErrorCheck( ::VecAssemblyEnd( vec ) ); }
    inline static void VecAXPY ( Vec y, PetscScalar alpha, Vec x) { ErrorCheck( ::VecAXPY( y, alpha, x ) ); }
    inline static void VecCopy ( Vec x, Vec y ) { ErrorCheck( ::VecCopy( x, y ) ); }

    template <class Comm>
    inline static void VecCreate ( const Comm& comm, Vec *vec )
    {
      ErrorCheck( ::VecCreate( castToPetscComm( comm ), vec ) );
    }

    template <class Comm>
    inline static void VecCreateGhost ( const Comm& comm, PetscInt n, PetscInt N, PetscInt nghost, const PetscInt ghosts[], Vec *vv )
      { ErrorCheck( ::VecCreateGhost( castToPetscComm( comm ), n, N, nghost, ghosts, vv ) ); }

    template <class Comm>
    inline static void VecCreateGhostBlock ( const Comm& comm, PetscInt bs, PetscInt n, PetscInt N, PetscInt nghost, const PetscInt ghosts[], Vec *vv )
    {
#if PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 2)
      std::unique_ptr< PetscInt[] > ghostsCopy( new PetscInt[ nghost * bs ] );
      for( PetscInt i = 0; i < nghost; ++i )
      {
        for( PetscInt j = 0; j < bs; ++j )
          ghostsCopy[ i*bs + j ] = ghosts[ i ]*bs + j;
      }
      VecCreateGhost( n, N, nghost * bs, ghostsCopy.get(), vv );
#else // #if PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 2)
      ErrorCheck( ::VecCreateGhostBlock( castToPetscComm( comm ), bs, n, N, nghost, ghosts, vv ) );
#endif // #else // #if PETSC_VERSION_MAJOR < 3 || (PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR <= 2)
    }

    inline static void VecDestroy ( Vec *v ) {
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::VecDestroy( *v ) );
      #else
        ErrorCheck( ::VecDestroy( v ) );
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    inline static void VecDot ( Vec x, Vec y, PetscScalar *val ) { ErrorCheck( ::VecDot( x, y, val ) ); }
    inline static void VecDuplicate ( Vec v, Vec *newv ) { ErrorCheck( ::VecDuplicate( v, newv ) ); }
    inline static void VecGetBlockSize ( Vec x, PetscInt* bs ) { ErrorCheck( ::VecGetBlockSize( x, bs ) ); }
    inline static void VecGetLocalSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetLocalSize( x, size ) ); }
    inline static void VecGetOwnershipRange ( Vec x, PetscInt *low, PetscInt *high ) { ErrorCheck( ::VecGetOwnershipRange( x, low, high ) ); }
    inline static void VecGetSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetSize( x, size ) ); }
    inline static void VecGetValues ( Vec x, PetscInt ni, const PetscInt ix[], PetscScalar y[] ) { ErrorCheck( ::VecGetValues( x, ni, ix, y ) ); }
    inline static void VecGhostGetLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostGetLocalForm( g, l ) ); }
    inline static void VecGhostRestoreLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostRestoreLocalForm( g, l ) ); }
    inline static void VecGhostUpdateBegin ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateBegin( g, insertmode, scattermode ) ); }
    inline static void VecGhostUpdateEnd ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateEnd( g, insertmode, scattermode ) ); }
    inline static void VecNorm ( Vec x, NormType type, PetscReal *val ) { ErrorCheck( ::VecNorm( x, type, val ) ); }
    inline static void VecScale ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecScale( x, alpha ) ); }
    inline static void VecSet ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecSet( x, alpha ) ); }
    inline static void VecSetBlockSize ( Vec x, PetscInt bs ) { ErrorCheck( ::VecSetBlockSize( x, bs ) ); }
    inline static void VecSetFromOptions ( Vec vec ) { ErrorCheck( ::VecSetFromOptions( vec ) ); }
    inline static void VecSetType ( Vec vec, const VecType method ) { ErrorCheck( ::VecSetType( vec, method ) ); }
    inline static void VecSetSizes ( Vec v, PetscInt n, PetscInt N ) { ErrorCheck( ::VecSetSizes( v, n, N ) ); }
    inline static void VecSetValue ( Vec v, int row, PetscScalar value, InsertMode mode ) { ErrorCheck( ::VecSetValue( v, row, value, mode ) ); }
    inline static void VecSetValuesBlocked ( Vec v, PetscInt ni, const PetscInt xi[], const PetscScalar values[], InsertMode mode ) { ErrorCheck( ::VecSetValuesBlocked( v, ni, xi, values, mode ) ); }
    inline static void VecView ( Vec vec, PetscViewer viewer ) { ErrorCheck( ::VecView( vec, viewer ) ); }

  } // namespace Petsc

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCCOMMON_HH
