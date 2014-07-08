# IGNORE_HEADERS
#CXXFLAGS=-g -Wall
# set default grid type for make all
GRIDTYPE=YASPGRID
#GRIDTYPE=PRISMGRID_ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_CUBE
#GRIDTYPE=UGGRID
GRIDDIM=3

EXTRA_DIST = checkleafcodim1.hh parameter *.dgf

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = cachingquadrature elementquadrature
TESTS = cachingquadrature elementquadrature

cachingquadrature_SOURCES = cachingquadrature.cc
elementquadrature_SOURCES = elementquadrature.cc

include $(top_srcdir)/am/global-rules

