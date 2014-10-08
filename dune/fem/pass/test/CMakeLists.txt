set(GRIDTYPE  YASPGRID)
set(GRIDDIM  2)
set(DIMRANGE  2)
set(POLORDER  2)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(SUBDIRS )
# test programs
set(TESTS  test-pass)
set(check_PROGRAMS  test-pass)
add_executable(test_pass test-pass.cc)
dune_target_link_libraries(test_pass ${DUNE_LIBS} ${LOCAL_LIBS})
add_dune_mpi_flags(test_pass )
set(DGFFILES  1dgrid.dgf 2dgrid.dgf 3dgrid.dgf)
set(EXTRA_DIST  ${DGFFILES} parameter)
set(CLEANFILES )
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
