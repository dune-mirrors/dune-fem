# IGNORE_HEADERS

EXTRA_DIST = cachequad_test.hh \
             albertagrid_fixture.hh \
             alugrid_fixture.hh \
             sgrid_fixture.hh \
             quad_test.hh 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DENABLE_PRIMSGRID -I$(HOME)/disk1/Dune/trunk/dune-prismgrid

check_PROGRAMS = tester cachingquadrature
tester_SOURCES = main.cc ../../misc/test.cc ../../misc/suite.cc \
	quad_test.cc cachequad_test.cc 

cachingquadrature_element_SOURCES = cachingquadrature.cc
cachingquadrature_element_CPPFLAGS = $(AM_CPPFLAGS)

include $(top_srcdir)/am/global-rules

