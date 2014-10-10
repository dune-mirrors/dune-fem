set(GRIDTYPE  YASPGRID)
set(GRIDDIM  2)
set(DIMRANGE  2)
set(POLORDER  2)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
add_definitions( "-DDIMRANGE=${DIMRANGE}" ) 
add_definitions( "-DPOLORDER=${POLORDER}" ) 

# test programs
set(TESTS  test_pass)
# copy data to build source to make tests work         

configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_executable(test_pass test-pass.cc)
dune_target_link_libraries(test_pass "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_pass )

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
