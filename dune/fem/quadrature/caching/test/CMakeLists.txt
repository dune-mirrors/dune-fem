# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
# also distribute the headers for the checks
set(EXTRA_DIST  cache_test.hh  refelem_test.hh  twist_test.hh  point_test.hh)
set(check_PROGRAMS  tester)
add_executable(tester main.cc ../../../misc/test.cc ../../../misc/suite.cc
        refelem_test.cc cache_test.cc twist_test.cc point_test.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
