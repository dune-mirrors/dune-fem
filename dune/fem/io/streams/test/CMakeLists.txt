set(TESTS  test_streams)

add_executable(test_streams test-streams.cc)
dune_target_link_libraries(test_streams "${DUNE_LIBS};${LOCAL_LIBS}")

set(streamstestdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/streams/test)
set(streamstest_HEADERS  test-streams.cc)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${streamstest_HEADERS} DESTINATION ${streamstestdir})

dune_add_test(${TESTS})
