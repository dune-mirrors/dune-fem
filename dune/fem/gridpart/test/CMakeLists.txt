if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

set( MARCOGRIDFILES 1dgrid.dgf  2dgrid.dgf 3dgrid.dgf )

exclude_from_headercheck( checkgridpart.hh )

# copy data to build source to make tests work
foreach( file ${MARCOGRIDFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME test_gplagrangeinterpolation SOURCES test-gplagrangeinterpolation.cc
COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM};COUNT_FLOPS;POLORDER=1"
LINK_LIBRARIES dunefem )

dune_add_test( NAME test_gridpart SOURCES test-gridpart.cc
  COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM};COUNT_FLOPS"
  LINK_LIBRARIES dunefem )
