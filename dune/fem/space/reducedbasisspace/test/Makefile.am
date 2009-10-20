LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)

POLORDER = 1

check_PROGRAMS = l2projection
l2projection_SOURCES = l2projection.cc

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -g -Wall -O0 -DPOLORDER=$(POLORDER)

testdir=$(includedir)/dune/fem/space/reducedbasisspace/test
test_HEADERS = l2projection.cc sinespace.hh
EXTRA_DIST = 2dgrid.dgf 3dgrid.dgf

include $(top_srcdir)/am/global-rules

