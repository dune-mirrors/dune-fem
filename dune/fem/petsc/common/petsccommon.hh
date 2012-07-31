// vim: set expandtab ts=2 sw=2 sts=2:
#ifndef DUNE_FEM_PETSCCOMMON_HH
#define DUNE_FEM_PETSCCOMMON_HH

/*
 * This file should be the first include in every dune-fem-petsc related header
 * It contains some common code in order to deal with PETSc in a consistent, object oriented
 * fashion
 */

#include <string>

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/mpimanager.hh>

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
    ErrorCode ErrorCheckHelper ( ErrorCode errorCode ) { CHKERRQ( errorCode ); return 0; }

    ErrorCode ErrorHandler ( MPI_Comm comm, int line, const char *function, const char *file, const char *dir, 
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

    void ErrorCheck ( ErrorCode errorCode )
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
    void initialize( int *argc, char ***args, const char file[], const char help[], bool ownHandler = true ) 
    {
      ::PetscInitialize( argc, args, file, help );

      if( ownHandler )
      {
        // set up our error handler
        if( Dune::MPIManager::rank() == 0 )
        {
          std::cerr << "INFORMATION: Setting up an own error handler for PETSc errors. If you want the default PETSc handler,\n"
                    << "INFORMATION: set the last argument of Dune::Petsc::initialize(...) to 'false'.\n";
        }
        ::PetscPushErrorHandler( &ErrorHandler, 0 );
      }
    }

    /*
     * This should be called just before the termination of the program.
     */
    void finalize ()
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
    void KSPCreate ( KSP *inksp ) { ErrorCheck( ::KSPCreate( FEM_PETSC_COMM_DEFAULT, inksp ) ); }
    void KSPDestroy ( KSP *ksp ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::KSPDestroy( *ksp ) ); 
      #else
        ErrorCheck( ::KSPDestroy( ksp ) ); 
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    void KSPGetPC ( KSP ksp, PC *pc ) { ErrorCheck( ::KSPGetPC( ksp, pc ) ); }
    void KSPSetFromOptions ( KSP ksp ) { ErrorCheck( ::KSPSetFromOptions( ksp ) ); }
    void KSPSetType ( KSP ksp, const KSPType type ) { ErrorCheck( ::KSPSetType( ksp, type ) ); }
    void KSPGMRESSetRestart ( KSP ksp, PetscInt restart ) { ErrorCheck( ::KSPGMRESSetRestart( ksp, restart ) ); }
    void KSPView ( KSP ksp, PetscViewer viewer = PETSC_VIEWER_STDOUT_(FEM_PETSC_COMM_DEFAULT)  ) 
    { 
      ErrorCheck( ::KSPView( ksp, viewer ) ); 
    }
    void KSPMonitorSet (KSP ksp, PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void*),
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
                               void *mctx,PetscErrorCode (*monitordestroy)(void*) 
      #else
                               void *mctx,PetscErrorCode (*monitordestroy)(void**)  
      #endif
                               )
    { 
        ErrorCheck( ::KSPMonitorSet( ksp, monitor, mctx, monitordestroy ) ); 
    }
    void KSPGetIterationNumber( KSP ksp, PetscInt* its ) 
    { ErrorCheck( ::KSPGetIterationNumber( ksp, its ) ); }

    void KSPSetOperators (KSP ksp, Mat Amat, Mat Pmat, MatStructure flag ) { ErrorCheck( ::KSPSetOperators( ksp, Amat, Pmat, flag ) ); }
    void KSPSetTolerances ( KSP ksp, PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits ) 
      { ErrorCheck( ::KSPSetTolerances( ksp, rtol, abstol, dtol, maxits ) ); }
    void KSPSolve ( KSP ksp, Vec b, Vec x ) { ErrorCheck( ::KSPSolve( ksp, b, x ) ); }
    void KSPSetPC ( KSP ksp, PC pc ) { ErrorCheck( ::KSPSetPC( ksp, pc ) ); }

    // preconditioning 
    void PCCreate  ( PC* pc) { ErrorCheck( ::PCCreate( FEM_PETSC_COMM_DEFAULT, pc ) ); }
    void PCDestroy ( PC* pc) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::PCDestroy( *pc ) ); 
      #else  
        ErrorCheck( ::PCDestroy( pc ) ); 
      #endif 
    }
    void PCSetType ( PC pc, const PCType type ) { ErrorCheck( ::PCSetType(  pc, type ) ); }
    void PCFactorSetLevels( PC pc, PetscInt level ) { ErrorCheck( ::PCFactorSetLevels(  pc, level ) ); }

    // matrix routines 
    void MatAssemblyBegin ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyBegin( mat, type ) ); }
    void MatAssemblyEnd ( Mat mat, MatAssemblyType type ) { ErrorCheck( ::MatAssemblyEnd( mat, type ) ); }
    void MatCreate ( Mat *A ) { ErrorCheck( ::MatCreate( FEM_PETSC_COMM_DEFAULT, A) ); }
    void MatCreateBlockMat ( Mat *A, PetscInt m, PetscInt n, PetscInt bs, PetscInt nz, PetscInt* nnz ) 
    { 
      ErrorCheck( ::MatCreateBlockMat( FEM_PETSC_COMM_DEFAULT, n, m, bs, nz, nnz, A) ); 
    }
    void MatDestroy ( Mat *A ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::MatDestroy( *A ) ); 
      #else
        ErrorCheck( ::MatDestroy( A ) ); 
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    }
    void MatGetOwnershipRange ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetOwnershipRange( mat, m, n ) ); }
    void MatGetSize ( Mat mat, PetscInt *m, PetscInt* n ) { ErrorCheck( ::MatGetSize( mat, m, n ) ); }
    void MatMult  ( Mat mat, Vec x, Vec y ) { ErrorCheck( ::MatMult( mat, x, y ) ); }
    void MatSetBlockSize ( Mat A, PetscInt bs ) { ErrorCheck( ::MatSetBlockSize( A, bs ) ); }
    void MatSetSizes ( Mat A, PetscInt m, PetscInt n, PetscInt M, PetscInt N ) { ErrorCheck( ::MatSetSizes( A, m, n, M, N ) ); }
    void MatSetFromOptions ( Mat B ) { ErrorCheck( ::MatSetFromOptions( B ) ); }
    void MatSetType ( Mat mat, const MatType matype ) { ErrorCheck( ::MatSetType( mat, matype ) ); }
    void MatSetValue ( Mat v, PetscInt i, PetscInt j, PetscScalar va, InsertMode mode ) { ErrorCheck( ::MatSetValue( v, i, j, va, mode ) ); }
    void MatSetValues ( Mat mat, PetscInt m, const PetscInt idxm[], PetscInt n, const PetscInt idxn[], const PetscScalar v[], InsertMode addv )
      { ErrorCheck( ::MatSetValues( mat, m, idxm, n, idxn, v, addv ) ); }
    void MatView ( Mat mat, PetscViewer viewer ) { ErrorCheck( ::MatView( mat, viewer ) ); }
    void MatZeroEntries ( Mat mat ) { ErrorCheck( ::MatZeroEntries( mat ) ); }
    void PetscBarrier ( PetscObject obj ) { ErrorCheck( ::PetscBarrier( obj ) ); }
    void PetscFinalize () { ErrorCheck( ::PetscFinalize() ); }
    void PetscInitialize( int *argc, char ***args, const char file[], const char help[] ) { ErrorCheck( ::PetscInitialize( argc, args, file, help ) ); }
    void PetscViewerASCIIOpen ( MPI_Comm comm, const char name[], PetscViewer *lab ) { ErrorCheck( ::PetscViewerASCIIOpen( comm, name, lab ) ); }
    void PetscViewerDestroy ( PetscViewer *viewer ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::PetscViewerDestroy( *viewer ) );  
      #else
        ErrorCheck( ::PetscViewerDestroy( viewer ) );  
      #endif
    }
    void PetscViewerSetFormat ( PetscViewer viewer, PetscViewerFormat format ) { ErrorCheck( ::PetscViewerSetFormat( viewer, format ) ); }
    void VecAssemblyBegin ( Vec vec ) { ErrorCheck( ::VecAssemblyBegin( vec ) ); }
    void VecAssemblyEnd ( Vec vec ) { ErrorCheck( ::VecAssemblyEnd( vec ) ); }
    void VecAXPY ( Vec y, PetscScalar alpha, Vec x) { ErrorCheck( ::VecAXPY( y, alpha, x ) ); }
    void VecCopy ( Vec x, Vec y ) { ErrorCheck( ::VecCopy( x, y ) ); }
    void VecCreate ( Vec *vec ) { ErrorCheck( ::VecCreate( FEM_PETSC_COMM_DEFAULT, vec ) ); }
    void VecCreateGhost ( PetscInt n, PetscInt N, PetscInt nghost, const PetscInt ghosts[], Vec *vv )
      { ErrorCheck( ::VecCreateGhost( FEM_PETSC_COMM_DEFAULT, n, N, nghost, ghosts, vv ) ); }
    void VecDestroy ( Vec *v ) { 
      #if PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
        ErrorCheck( ::VecDestroy( *v ) ); 
      #else
        ErrorCheck( ::VecDestroy( v ) ); 
      #endif // PETSC_PETSC_VERSION_MAJOR <= 3 && PETSC_VERSION_MINOR < 2
    } 
    void VecDot ( Vec x, Vec y, PetscScalar *val ) { ErrorCheck( ::VecDot( x, y, val ) ); }
    void VecDuplicate ( Vec v, Vec *newv ) { ErrorCheck( ::VecDuplicate( v, newv ) ); }
    void VecGetLocalSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetLocalSize( x, size ) ); }
    void VecGetOwnershipRange ( Vec x, PetscInt *low, PetscInt *high ) { ErrorCheck( ::VecGetOwnershipRange( x, low, high ) ); }
    void VecGetSize ( Vec x, PetscInt *size ) { ErrorCheck( ::VecGetSize( x, size ) ); }
    void VecGetValues ( Vec x, PetscInt ni, const PetscInt ix[], PetscScalar y[] ) { ErrorCheck( ::VecGetValues( x, ni, ix, y ) ); }
    void VecGhostGetLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostGetLocalForm( g, l ) ); }
    void VecGhostRestoreLocalForm ( Vec g, Vec *l ) { ErrorCheck( ::VecGhostRestoreLocalForm( g, l ) ); }
    void VecGhostUpdateBegin ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateBegin( g, insertmode, scattermode ) ); }
    void VecGhostUpdateEnd ( Vec g, InsertMode insertmode, ScatterMode scattermode ) { ErrorCheck( ::VecGhostUpdateEnd( g, insertmode, scattermode ) ); }
    void VecNorm ( Vec x, NormType type, PetscReal *val ) { ErrorCheck( ::VecNorm( x, type, val ) ); }
    void VecScale ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecScale( x, alpha ) ); }
    void VecSet ( Vec x, PetscScalar alpha ) { ErrorCheck( ::VecSet( x, alpha ) ); }
    void VecSetFromOptions ( Vec vec ) { ErrorCheck( ::VecSetFromOptions( vec ) ); }
    void VecSetType ( Vec vec, const VecType method ) { ErrorCheck( ::VecSetType( vec, method ) ); }
    void VecSetSizes ( Vec v, PetscInt n, PetscInt N ) { ErrorCheck( ::VecSetSizes( v, n, N ) ); }
    void VecSetValue ( Vec v, int row, PetscScalar value, InsertMode mode ) { ErrorCheck( ::VecSetValue( v, row, value, mode ) ); }
    void VecView ( Vec vec, PetscViewer viewer ) { ErrorCheck( ::VecView( vec, viewer ) ); }
      
  } // namespace Petsc

} // namespace Dune

#endif // #if HAVE_PETSC

#endif // DUNE_FEM_PETSCCOMMON_HH
