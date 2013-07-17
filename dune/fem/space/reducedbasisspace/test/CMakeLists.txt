set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)
set(POLORDER  1)
# removed until somebody who actually needs this fixes it
#check_PROGRAMS = l2projection
add_executable(l2projection l2projection.cc)
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -g -Wall -O0 -DPOLORDER=${POLORDER})
set(testdir ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/reducedbasisspace/test)
set(test_HEADERS  l2projection.cc sinespace.hh)
set(EXTRA_DIST  2dgrid.dgf 3dgrid.dgf)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${test_HEADERS} DESTINATION ${testdir})