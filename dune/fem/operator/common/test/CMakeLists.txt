add_executable(mappings mappingtest.cc)
dune_target_link_libraries(mappings "${DUNE_LIBS};${LOCAL_LIBS}")

dune_fem_add_test(mappings)
