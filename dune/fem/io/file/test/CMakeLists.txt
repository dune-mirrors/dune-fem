set(POLORDER  2)
set(WANT_GRAPE  0)
# set default grid type for make all
set(GRIDTYPE YASPGRID)
set(GRIDDIM 2)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS} ${MPI_LDFLAGS})
# AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
set(check_PROGRAMS  persistencemanagertest dataoutputtest)
add_executable(persistencemanagertest persistencemanagertest.cc)
add_executable(dataoutputtest dataoutputtest.cc)
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -DPOLORDER=${POLORDER} -DWANT_GRAPE=${WANT_GRAPE})
set(CXXFLAGS  -g)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
