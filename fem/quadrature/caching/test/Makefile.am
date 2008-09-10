LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

cachingtestdir =  $(includedir)/dune/fem/quadrature/caching/test
cachingtest_HEADERS = albertagrid_fixture.hh  cache_test.hh \
refelem_test.hh   twist_test.hh alugrid_fixture.hh \
point_test.hh  sgrid_fixture.hh

check_PROGRAMS = tester 
tester_SOURCES = main.cc ../../../misc/test.cc ../../../misc/suite.cc \
	refelem_test.cc cache_test.cc twist_test.cc point_test.cc

include $(top_srcdir)/am/global-rules

