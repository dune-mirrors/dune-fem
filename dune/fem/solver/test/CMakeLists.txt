# copy data to build source to make tests work
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)

add_executable( odetest odetest.cc)
dune_target_link_libraries( odetest "${DUNE_LIBS};${LOCAL_LIBS}")

add_executable( test_timeprovider test-timeprovider.cc)
dune_target_link_libraries( test_timeprovider "${DUNE_LIBS};${LOCAL_LIBS}")

dune_fem_add_test( odetest test_timeprovider )
