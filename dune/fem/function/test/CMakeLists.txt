set( DATAFILES
  1dgrid.dgf
  2dgrid.dgf
  3dgrid.dgf
)

foreach( datafile ${DATAFILES} )
  configure_file( ${datafile} ${CMAKE_CURRENT_BINARY_DIR}/${datafile} COPYONLY )
endforeach()

include( GridSelector )
add_dune_gridselector( gridselector.hh )

add_executable( dftest EXCLUDE_FROM_ALL dftest.cc )
target_depend_dune_gridselector( dftest )
dune_target_link_libraries( dftest "${DUNE_LIBS};${LOCAL_LIBS}" )

dune_add_test( dftest )
