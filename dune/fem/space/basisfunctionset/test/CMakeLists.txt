if( NOT GRIDTYPE )
  set(GRIDTYPE  YASPGRID)
endif()
if( NOT GRIDDIM )
  set(GRIDDIM  2)
endif()

# copy data to build source to make tests work
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)

set( BASISFUNCTIONSETS vectorialbasisfunctionset defaultbasisfunctionset  simplebasisfunctionset tuplebasisfunctionset )
foreach( bset ${BASISFUNCTIONSETS} )
  dune_add_test( NAME test_${bset} SOURCES test-${bset}.cc
  COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM};DIMRANGE=5;POLORDER=3"
  LINK_LIBRARIES dunefem )
endforeach()

set_property(TARGET test_vectorialbasisfunctionset APPEND PROPERTY COMPILE_DEFINITIONS "USE_VERTICAL_DOF_ALIGNMENT=1" )

dune_install(checkbasisfunctionset.hh)
