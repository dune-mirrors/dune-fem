if( NOT GRIDTYPE )
  set(GRIDTYPE YASPGRID)
endif()
if( NOT GRIDDIM )
  set(GRIDDIM 2)
endif()

# copy data to build source to make tests work
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)


set( DEFAULTFLAGS "-D${GRIDTYPE}" "-DGRIDDIM=${GRIDDIM}" "-DDIMENSION=2" "-DPOLORDER=2" "-DTOPOLOGYTYPE=Dune::GeometryTypes::simplex" "-DSCALAR=0" )

if( ${FEM_TORTURE_TESTS} )
  # dune_add_test( NAME tester SOURCES main.cc
  # COMPILE_DEFINITIONS  ${DEFAULTFLAGS}
  # LINK_LIBRARIES dunefem )

  dune_add_test( NAME generic SOURCES generic.cc
  COMPILE_DEFINITIONS  ${DEFAULTFLAGS}
  LINK_LIBRARIES dunefem )
endif()
