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

dune_add_test( NAME dftest SOURCES dftest.cc
COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}"
LINK_LIBRARIES dunefem )
