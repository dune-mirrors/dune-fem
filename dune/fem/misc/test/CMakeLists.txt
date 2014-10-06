set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(SUBDIRS )

set(EXTRA_DIST )
set(CLEANFILES )

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
#add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})
