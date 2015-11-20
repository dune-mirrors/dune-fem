if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 3 )
endif()

# copy data to build source to make tests work
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

# the default compiler flags
set( DEFAULTFLAGS "GRIDDIM=${GRIDDIM};${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0" )

set( TESTS
  adapt
  arraytest
  combinedspacetest
  functiontupleadapt
  hierarchicspace
  interpolation
  l2projection
  lagrangeadapt
  lagrangeglobalrefine
  lagrangeinterpolation
  legendretest
  padapt
)

# create test targets, for standart tests
foreach( test ${TESTS} )
  dune_add_test( NAME ${test} SOURCES ${test}.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}"
  LINK_LIBRARIES dunefem )
endforeach()


# dgcomm tests should run on GRIDDIM=3 meshes
dune_add_test( NAME dgcomm SOURCES dgcomm.cc
COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0;WANT_CACHED_COMM_MANAGER=0"
LINK_LIBRARIES dunefem MPI_RANKS 1 3 4 TIMEOUT 9999999 )

dune_add_test( NAME dgcomm_cached SOURCES dgcomm.cc
COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0"
LINK_LIBRARIES dunefem MPI_RANKS 1 3 4 TIMEOUT 9999999 )


# test adaptation of combinedspaces
dune_add_test( NAME adapt_combinedspace_pointbased SOURCES adapt.cc
COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_COMBINED_SPACE;POINTBASED"
LINK_LIBRARIES dunefem )

dune_add_test( NAME adapt_combinedspace_variablebased SOURCES adapt.cc
COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_COMBINED_SPACE"
LINK_LIBRARIES dunefem )
