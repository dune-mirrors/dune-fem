set( POLORD  3 )
set( GRIDTYPE  YASPGRID )
set( GRIDDIM  2 )

# remove this until this test is fixed again
#set( ADAPTPROG  prod )
# remove if-condition after test is fixed again
if( ${ADAPTPROG} )
  add_definitions( "-D${GRIDTYPE}" )
  add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
  add_definitions( "-DPOLORD=${POLORD}" )
  
  # copy data to build source to make tests work         
  configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)

  add_executable(prod productfunctiontest.cc)
  dune_target_link_libraries(prod "${DUNE_LIBS};${LOCAL_LIBS}")
  
  dune_add_test(prod)
endif( ${ADAPTPROG} )
