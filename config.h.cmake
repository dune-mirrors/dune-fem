/* begin dune-fem */

/* Define to the version of dune-fem */
#define DUNE_FEM_VERSION "${DUNE_FEM_VERSION}"

/* Define to the major version of dune-fem */
#define DUNE_FEM_VERSION_MAJOR ${DUNE_FEM_VERSION_MAJOR}

/* Define to the minor version of dune-fem */
#define DUNE_FEM_VERSION_MINOR ${DUNE_FEM_VERSION_MINOR}

/* Define to the revision of dune-fem */
#define DUNE_FEM_VERSION_REVISION ${DUNE_FEM_VERSION_REVISION}


/* Define if we have name for endian header (including brackets) */
#cmakedefine SYSTEM_ENDIAN_HEADER "${SYSTEM_ENDIAN_HEADER}"

/* Define actual name of xdr_uint_64_t method */
#cmakedefine XDR_UINT64_FUNC ${XDR_UINT64_FUNC}

/* Define if we have pthreads */
#cmakedefine HAVE_PTHREAD 1

/* Define if we have papi */
#cmakedefine HAVE_PAPI 1

/* Define if we have eigen */
#cmakedefine HAVE_EIGEN 1

/* Define if we have sionlib */
#cmakedefine HAVE_SIONLIB 1

/* Define if we have PETSc */
#cmakedefine HAVE_PETSC ENABLE_PETSC

/* Define if we have AmgXSolver */
#cmakedefine HAVE_AMGXSOLVER 1

/* Define if we have ViennaCL */
#cmakedefine HAVE_VIENNACL 1

/* Define if we have OpenCL */
#cmakedefine HAVE_OPENCL 1

/* end dune-fem */
