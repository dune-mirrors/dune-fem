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

dune_add_test( NAME lagrangeinterpolation_dynamic SOURCES lagrangeinterpolation.cc
               COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )
set( TESTS ${TESTS} lagrangeinterpolation_dynamic )

# dgl2projection tests for different discrete functions and spaces
dune_add_test( NAME dgl2projection SOURCES dgl2projection.cc
COMPILE_DEFINITIONS "POLORDER=2;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} dgl2projection )

dune_add_test( NAME vtxprojection_adaptive SOURCES vtxprojection.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} vtxprojection_adaptive )

dune_add_test( NAME adaptation SOURCES adaptation.cc
  COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;ALUGRID_CONFORM" MPI_RANKS 1 2 3 4 TIMEOUT 9999999 CMAKE_GUARD dune-alugrid_FOUND)
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

dune_add_test( NAME l2projection_adaptive SOURCES l2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} l2projection_adaptive )

dune_add_test( NAME l2projection_istl SOURCES l2projection.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} l2projection_istl )

if( ${FEM_TORTURE_TESTS} )
  dune_add_test( NAME l2projection_istl_complex SOURCES l2projection.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;USE_COMPLEX;${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} l2projection_istl_complex )

  if( ${SuiteSparse_UMFPACK_FOUND} )
    dune_add_test( NAME l2projection_adaptive_umfpack SOURCES l2projection.cc
    COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS};USE_UMFPACK" )
    set( TESTS ${TESTS} l2projection_adaptive_umfpack )
  endif()

  dune_add_test( NAME test-tupleoperator SOURCES test-tupleoperator.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;USE_COMPLEX;${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} test-tupleoperator )
endif()

dune_add_test( NAME test-periodic SOURCES test-periodic.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}" CMAKE_GUARD dune-alugrid_FOUND)
set( TESTS ${TESTS} test-periodic )

dune_add_test( NAME referencesolution SOURCES referencesolution.cc
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
