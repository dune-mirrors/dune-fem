# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

EXTRA_DIST = cachequad_test.hh \
             albertagrid_fixture.hh \
             alugrid_fixture.hh \
             sgrid_fixture.hh \
             quad_test.hh 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = tester cachingquadrature

tester_SOURCES = main.cc ../../misc/test.cc ../../misc/suite.cc quad_test.cc cachequad_test.cc

cachingquadrature_SOURCES = cachingquadrature.cc

include $(top_srcdir)/am/global-rules

