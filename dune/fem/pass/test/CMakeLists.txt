add_executable(test_insertoperatorpass test-insertoperatorpass.cc)
dune_target_link_libraries(test_insertoperatorpass "${DUNE_LIBS};${LOCAL_LIBS}")
add_dune_mpi_flags(test_insertoperatorpass)
dune_add_test(test_insertoperatorpass)
