if( NOT GRIDTYPE )
  set( GRIDTYPE YASPGRID )
endif()

if( NOT GRIDDIM )
  set( GRIDDIM 2 )
endif()

configure_file( 2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY )

dune_add_test( NAME mappingtest SOURCES mappingtest.cc LINK_LIBRARIES dunefem )
dune_add_test( NAME dofconstraintest SOURCES dofconstraintest.cc LINK_LIBRARIES dunefem
               COMPILE_DEFINITIONS "${GRIDTYPE}" "GRIDDIM=${GRIDDIM}")
