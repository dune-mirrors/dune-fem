# if no grid type is specified during configuration time we choose dune-alugrid simplex
if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

# same as for GridType, we set the default grid dimension to 2
if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

# copy data to build source to make tests work
set( TESTDATAFILES 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf parameter )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

# some default flags needed very often
set( DEFAULTFLAGS "COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

# variable which gathers all test targets
set( TESTS )

# lagrangeinterpolation tests
if( ${GRIDDIM} LESS 3 )
  set( MAXORDER 4 )
else()
  set( MAXORDER 2 )
endif()

foreach( order RANGE 1 ${MAXORDER} )
  dune_add_test( NAME lagrangeinterpolation_p${order} SOURCES lagrangeinterpolation.cc
                 COMPILE_DEFINITIONS "POLORDER=${order};${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} lagrangeinterpolation_p${order} )
endforeach()
dune_add_test( NAME lagrangeinterpolation_dynamic SOURCES lagrangeinterpolation.cc
               COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )
set( TESTS ${TESTS} lagrangeinterpolation_dynamic )

# dgl2projection tests for different discrete functions
dune_add_test( NAME dgl2projection_adaptive SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection_adaptive )

set( DFS blockvectorfunction vectorfunction blockvectordiscretefunction )
foreach( df ${DFS} )
  string( TOUPPER ${df} DF )
  dune_add_test( NAME dgl2projection_${df} SOURCES dgl2projection.cc
                 COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_${DF};${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} dgl2projection_${df} )
endforeach()

dune_add_test( NAME dgl2projection_vector_int SOURCES dgl2projection.cc
               COMPILE_DEFINITIONS "POLORDER=1;USE_VECTORFUNCTION;USE_DOFTYPE_INT;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection_vector_int )

# only run this test for YaspGrid
dune_add_test( NAME dgl2projection_legendre SOURCES dgl2projection.cc
               COMPILE_DEFINITIONS "POLORDER=1;YASPGRID;USE_LEGENDRESPACE;GRIDDIM=${GRIDDIM}" )
set( TESTS ${TESTS} dgl2projection_legendre )

dune_add_test( NAME dgl2projection_hierarchiclegendre SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_HIERARCHICLEGENDRESPACE;YASPGRID;GRIDDIM=${GRIDDIM}" )
set( TESTS ${TESTS} dgl2projection_hierarchiclegendre )

dune_add_test( NAME dgl2projection_fv SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "USE_FVSPACE;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection_fv )


dune_add_test( NAME dgl2projection_filtered SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_FILTEREDGRID;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection_filtered )

dune_add_test( NAME dgl2projection_idgridpart SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_IDGRIDPART;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection_idgridpart )


dune_add_test( NAME vtxprojection_adaptive SOURCES vtxprojection.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} vtxprojection_adaptive )

dune_add_test( NAME combinedspaceadaptation SOURCES adaptation.cc
COMPILE_DEFINITIONS "COMBINEDSPACE;GRIDDIM=3;WORLDDIM=3;${GRIDTYPE}" )
set( TESTS ${TESTS} combinedspaceadaptation )

dune_add_test( NAME adaptation SOURCES adaptation.cc
COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE}" MPI_RANKS 1 2 3 4 TIMEOUT 9999999 )
set( TESTS ${TESTS} adaptation )

dune_add_test(
  NAME l2projection_petsc
  SOURCES l2projection.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_PETSCDISCRETEFUNCTION;${DEFAULTFLAGS}"
  LINK_LIBRARIES dunefem
  MPI_RANKS 1 2 3 4
  TIMEOUT 300
  CMAKE_GUARD HAVE_PETSC
)

if( ${SuiteSparse_UMFPACK_FOUND} )
  dune_add_test( NAME l2projection_adaptive_umfpack SOURCES l2projection.cc
  COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS};USE_UMFPACK" )
  set( TESTS ${TESTS} l2projection_adaptive_umfpack )
endif()

dune_add_test( NAME l2projection_adaptive SOURCES l2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} l2projection_adaptive )

dune_add_test( NAME l2projection_istl SOURCES l2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} l2projection_istl )

dune_add_test( NAME l2projection_istl_complex SOURCES l2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;USE_COMPLEX;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} l2projection_istl_complex )


dune_add_test( NAME test-tupleoperator SOURCES test-tupleoperator.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;USE_COMPLEX;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} test-tupleoperator )

if( HAVE_EIGEN )
  dune_add_test( NAME l2projection_eigen SOURCES l2projection.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_EIGEN;${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} l2projection_eigen )
endif()

dune_add_test( NAME referencesolution SOURCES referencesolution
COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_HIERARCHICLEGENDRESPACE;YASPGRID;GRIDDIM=${GRIDDIM}" )
set( TESTS ${TESTS} referencesolution )

#finally add all pkg flags and local libs to tests
foreach( test ${TESTS} )
  target_link_dune_default_libraries( ${test} )
  dune_target_link_libraries( ${test} dunefem )
endforeach()

exclude_from_headercheck( dfspace.hh testgrid.hh )

dune_add_test(SOURCES test-intersectionindexset.cc LINK_LIBRARIES dunefem CMAKE_GUARD dune-alugrid_FOUND)

dune_add_test(
  SOURCES test-hierarchicall2projection.cc
  COMPILE_DEFINITIONS "YASPGRID" "GRIDDIM=2"
  LINK_LIBRARIES dunefem
  CMAKE_GUARD dune-istl_FOUND
)
