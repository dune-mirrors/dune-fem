# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)
set(EXTRA_DIST  testgrid.hh)
set(AM_CXXFLAGS -g -Wall)
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${MPI_LDFLAGS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  -Wall ${ALL_PKG_CPPFLAGS})
set(check_PROGRAMS  test)
add_executable(test test.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
