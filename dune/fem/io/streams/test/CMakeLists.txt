# IGNORE_HEADERS
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(AM_CXXFLAGS  -O0 -g -Wall ${ALL_PKG_CXXFLAGS})
set(TESTS  test-streams)
set(check_PROGRAMS  test-streams)
add_executable(test_streams test-streams.cc)
dune_target_link_libraries(test_streams ${DUNE_LIBS})
set(streamstestdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/streams/test)
set(streamstest_HEADERS  test-streams.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${streamstest_HEADERS} DESTINATION ${streamstestdir})
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
