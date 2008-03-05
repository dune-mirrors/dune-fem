LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

DIMENSION = 2
POLORDER = 2
GEOMETRYTYPE = simplex

check_PROGRAMS = tester generic 
tester_SOURCES = main.cc ../../../misc/test.cc ../../../misc/suite.cc 
generic_SOURCES = generic.cc
generic_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMENSION=$(DIMENSION) -DPOLORDER=$(POLORDER) -DGEOMETRYTYPE=$(GEOMETRYTYPE)

testdir=$(includedir)/dune/fem/space/lagrangespace/test
test_HEADERS = basefunctiontest.hh  mappertest.hh basefunctiontest.cc \
mappertest.cc
EXTRA_DIST = 2dgrid.dgf cube.dgf  onesimplex.dgf
