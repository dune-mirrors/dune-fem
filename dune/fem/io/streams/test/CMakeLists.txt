add_executable(test_streams test-streams.cc)
dune_target_link_libraries(test_streams "${DUNE_LIBS};${LOCAL_LIBS}")

dune_install(test-streams.cc)

dune_fem_add_test(test_streams)
