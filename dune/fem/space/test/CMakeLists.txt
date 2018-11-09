if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

# copy data to build source to make tests work
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(2dgrid_8.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid_8.dgf COPYONLY)
configure_file(3dgrid_8.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid_8.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

# the default compiler flags
set( DEFAULTFLAGS "GRIDDIM=${GRIDDIM};${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0" )

set( TESTS
  adapt
  lagrangeadapt
  padapt
)

set( EXPENSIVE_TESTS
  combinedspacetest
  functiontupleadapt
  hierarchicspace
  interpolation
  l2projection
  lagrangeglobalrefine
  lagrangeinterpolation
  legendretest
  localadapter
)

if( ${TORTURE_TESTS} )
  foreach(i ${EXPENSIVE_TESTS})
    list(APPEND TESTS ${i})
  endforeach()
endif()

# create test targets, for standart tests
foreach( test ${TESTS} )
  dune_add_test( NAME ${test} SOURCES ${test}.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}"
  LINK_LIBRARIES dunefem )
endforeach()


# dgcomm tests should run on GRIDDIM=3 meshes
dune_add_test( NAME dgcomm SOURCES dgcomm.cc
COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0;WANT_CACHED_COMM_MANAGER=0"
LINK_LIBRARIES dunefem MPI_RANKS 1 2 3 4 TIMEOUT 9999999 )

dune_add_test( NAME dgcomm_cached SOURCES dgcomm.cc
COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE};DIMRANGE=3;POLORDER=2;WANT_GRAPE=0"
LINK_LIBRARIES dunefem MPI_RANKS 1 2 3 4 TIMEOUT 9999999 )

dune_add_test( SOURCES test-slavedofs.cc LINK_LIBRARIES dunefem MPI_RANKS 1 2 4 TIMEOUT 300 COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_COMBINED_SPACE" )
dune_add_test( SOURCES test-raviartthomasinterpolation.cc LINK_LIBRARIES dunefem COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )

if( ${TORTURE_TESTS} )
  # test adaptation of combinedspaces
  dune_add_test( NAME adapt_combinedspace_pointbased SOURCES adapt.cc
  COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_COMBINED_SPACE;POINTBASED"
  LINK_LIBRARIES dunefem )

  dune_add_test( NAME adapt_combinedspace_variablebased SOURCES adapt.cc
  COMPILE_DEFINITIONS "${DEFAULTFLAGS};USE_COMBINED_SPACE"
  LINK_LIBRARIES dunefem )

  dune_add_test( SOURCES test-localadaptation.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}" LINK_LIBRARIES dunefem MPI_RANKS 1 2 3 4 TIMEOUT 9999999 )
  dune_add_test( SOURCES test-globalrefine.cc COMPILE_DEFINITIONS "${DEFAULTFLAGS}" LINK_LIBRARIES dunefem )

  dune_add_test( SOURCES test-balladapt.cc LINK_LIBRARIES dunefem MPI_RANKS 1 2 4 TIMEOUT 300 )

  dune_add_test( SOURCES test-tupledgspace.cc LINK_LIBRARIES dunefem COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )


  dune_add_test( SOURCES test-spaceinterpolation.cc LINK_LIBRARIES dunefem COMPILE_DEFINITIONS "${DEFAULTFLAGS}" )
endif()

install( FILES checklocalinterpolation.hh DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/test )
