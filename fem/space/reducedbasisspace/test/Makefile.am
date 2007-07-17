LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -g -Wall -O0

check_PROGRAMS = l2projection
l2projection_SOURCES = l2projection.cc

testdir=$(includedir)/fem/space/reducedbasisspace/test
test_HEADERS = l2projection.cc
EXTRA_DIST = 2dgrid.dgf 3dgrid.dgf
