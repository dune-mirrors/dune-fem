set( GRIDTYPE YASPGRID )
set( GRIDDIM 2 )

set( DATAFILES
  1dgrid.dgf
  2dgrid.dgf
  3dgrid.dgf
)

foreach( datafile ${DATAFILES} )
  configure_file( ${datafile} ${CMAKE_CURRENT_BINARY_DIR}/${datafile} COPYONLY )
endforeach()

add_executable( dftest EXCLUDE_FROM_ALL dftest.cc )
dune_target_link_libraries( dftest "${DUNE_LIBS};${LOCAL_LIBS}" )
target_compile_definitions( dftest PRIVATE "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}" )

dune_fem_add_test( dftest )
