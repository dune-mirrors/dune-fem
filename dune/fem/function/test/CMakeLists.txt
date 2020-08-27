if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

set( DATAFILES
  1dgrid.dgf
  2dgrid.dgf
  3dgrid.dgf
)

foreach( datafile ${DATAFILES} )
  configure_file( ${datafile} ${CMAKE_CURRENT_BINARY_DIR}/${datafile} COPYONLY )
endforeach()

dune_add_test( NAME test_filtered_df SOURCES test-filtered-df.cc
COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME test_tuplediscretefunction SOURCES test-tuplediscretefunction.cc
COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME test_df SOURCES test-df.cc
COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test( NAME test_lfa SOURCES test-lfa.cc
COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )

dune_add_test(
  SOURCES test-hierarchicalfunction.cc
  CMAKE_GUARD dune-istl_FOUND
  COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
  LINK_LIBRARIES dunefem
)

dune_add_test(
  SOURCES test-hierarchicalfunction-comm.cc
  CMAKE_GUARD dune-istl_FOUND
  COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
  LINK_LIBRARIES dunefem
  MPI_RANKS 2 3 4
  TIMEOUT 300
)
