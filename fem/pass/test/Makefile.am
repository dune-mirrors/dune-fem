# $Id$

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 
#AM_LDFLAGS = $(top_builddir)/lib/libdune.la

check_PROGRAMS = tester 
tester_SOURCES = main.cc test.cc suite.cc pass_test.cc helper_test.cc \
	test.hh suite.hh pass_test.hh helper_test.hh passstub.hh dgstub.hh
#tester_CXXFLAGS = $(ALUGRID_CPPFLAGS) $(AM_CPPFLAGS)
#tester_LDADD = $(ALUGRID_LDFLAGS) $(ALUGRID_LIBS) $(AM_LDFLAGS)

