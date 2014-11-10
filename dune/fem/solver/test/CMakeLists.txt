add_executable(odetest odetest.cc)
dune_target_link_libraries(odetest "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(odesystemtest odesystemtest.cc)
dune_target_link_libraries(odesystemtest "${DUNE_LIBS};${LOCAL_LIBS}")

dune_add_test(odetest odesystemtest)
