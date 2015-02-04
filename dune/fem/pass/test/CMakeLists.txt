set(GRIDTYPE  YASPGRID)
set(GRIDDIM  2)
set(DIMRANGE  2)
set(POLORDER  2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DDIMRANGE=${DIMRANGE}" ) 
add_definitions( "-DPOLORDER=${POLORDER}" ) 

# copy data to build source to make tests work         
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_executable(test_pass test-pass.cc)
dune_target_link_libraries(test_pass "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_pass )
dune_add_test(test_pass)

add_executable(test_insertoperatorpass test-insertoperatorpass.cc)
dune_target_link_libraries(test_insertoperatorpass "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_insertoperatorpass )
dune_add_test(test_insertoperatorpass)
