# if no grid type is specified during configuration time we choose dune-alugrid simplex
if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

# same as for GridType, we set the default grid dimension to 2
if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()
set( DEFAULTFLAGS "COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

# copy data to build source to make tests work
set( TESTDATAFILES 2dgrid.dgf )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME unitrows SOURCES unitrows.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;${DEFAULTFLAGS}"
  LINK_LIBRARIES dunefem )
set( TESTS ${TESTS} unitrows )

if( HAVE_PETSC )
  dune_add_test( NAME laplace_petsc SOURCES laplace.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_PETSCDISCRETEFUNCTION;${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} laplace_petsc )
endif()

dune_add_test( NAME laplace_adaptive SOURCES laplace.cc
COMPILE_DEFINITIONS "POLORDER=1;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} laplace_adaptive )

dune_add_test( NAME laplace_istl SOURCES laplace.cc
COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;${DEFAULTFLAGS}" )
set( TESTS ${TESTS} laplace_istl )

# dune_add_test( NAME test-tupleoperator SOURCES test-tupleoperator.cc
# COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORFUNCTION;USE_COMPLEX;${DEFAULTFLAGS}" )
# set( TESTS ${TESTS} test-tupleoperator )

if( HAVE_EIGEN )
  dune_add_test( NAME laplace_eigen SOURCES laplace.cc
  COMPILE_DEFINITIONS "POLORDER=1;USE_EIGEN;${DEFAULTFLAGS}" )
  set( TESTS ${TESTS} laplace_eigen )
endif()

#finally add all pkg flags and local libs to tests
foreach( test ${TESTS} )
  target_link_dune_default_libraries( ${test} )
  dune_target_link_libraries( ${test} dunefem )
endforeach()
