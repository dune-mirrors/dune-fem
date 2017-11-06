# if no grid type is specified during configuration time we choose yaspgrid
if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

# same as for GridType, we set the default grid dimension to 2
if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

# copy data to build source to make tests work
set( TESTDATAFILES 2dgrid.dgf 3dgrid.dgf parameter )
foreach( file ${TESTDATAFILES} )
  configure_file( ${file} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY )
endforeach()

dune_add_test( NAME backuprestore SOURCES backuprestore.cc
  COMPILE_DEFINITIONS "POLORDER=2;WANT_GRAPE=0;${GRIDTYPE};GRIDDIM=${GRIDDIM};DATA_PATH=\"${CMAKE_CURRENT_BINARY_DIR}\"" LINK_LIBRARIES dunefem MPI_RANKS 1 2 3 4 TIMEOUT 9999999 )

dune_add_test( NAME persistencemanagertest SOURCES persistencemanagertest.cc
COMPILE_DEFINITIONS "POLORDER=2;WANT_GRAPE=0;${GRIDTYPE};GRIDDIM=${GRIDDIM}" LINK_LIBRARIES dunefem )

dune_add_test( NAME dataoutputtest SOURCES dataoutputtest.cc
COMPILE_DEFINITIONS "POLORDER=2;WANT_GRAPE=0;${GRIDTYPE};GRIDDIM=${GRIDDIM}" LINK_LIBRARIES dunefem )
