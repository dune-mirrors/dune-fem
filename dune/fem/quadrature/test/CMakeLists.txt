# IGNORE_HEADERS
#CXXFLAGS=-g -Wall
# set default grid type for make all
set(GRIDTYPE YASPGRID)
#GRIDTYPE=PRISMGRID_ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_CUBE
#GRIDTYPE=UGGRID
set(GRIDDIM 3)
set(EXTRA_DIST  checkleafcodim1.hh)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(check_PROGRAMS  cachingquadrature)
add_executable(cachingquadrature cachingquadrature.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
