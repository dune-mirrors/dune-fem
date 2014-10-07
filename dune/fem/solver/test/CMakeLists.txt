set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})

set(TESTS  odetest odesystemtest)
set(check_PROGRAMS  odetest odesystemtest)

add_executable(odetest odetest.cc)
dune_target_link_libraries(odetest ${DUNE_LIBS} )
add_executable(odesystemtest odesystemtest.cc)
dune_target_link_libraries(odesystemtest ${DUNE_LIBS} )

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
