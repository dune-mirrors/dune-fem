if( NOT GRIDTYPE )
  set(GRIDTYPE SGRID)
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

set( MARCOGRIDFILES 1dgrid.dgf  2dgrid.dgf 3dgrid.dgf )

# copy data to build source to make tests work
foreach( file ${MARCOGRIDFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME gp_lagrangeinterpolation SOURCES lagrangeinterpolation.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM};COUNT_FLOPS;POLORDER=1"
LINK_LIBRARIES dunefem )

set( gp filteredgridpart idgridpart geogridpart )
foreach( gp ${GRIDPARTS} )
  dune_add_test( NAME test_${gp} SOURCES test-${GP}.cc
  COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM};COUNT_FLOPS"
  LINK_LIBRARIES dunefem )
endforeach()
