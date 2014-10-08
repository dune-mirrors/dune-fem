LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

DIMENSION = 2
POLORDER = 2
GEOMETRYTYPE = simplex
SCALAR = 0

check_PROGRAMS = tester generic testnonscalarmapper 
TESTS = tester generic testnonscalarmapper 

tester_SOURCES = main.cc 
generic_SOURCES = generic.cc
generic_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMENSION=$(DIMENSION) -DPOLORDER=$(POLORDER) -DGEOMETRYTYPE=$(GEOMETRYTYPE)

testnonscalarmapper_SOURCES = testnonscalarmapper.cc
testnonscalarmapper_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DSCALAR=$(SCALAR)


testdir=$(includedir)/dune/fem/space/lagrange/test
test_HEADERS = basisfunctiontest.hh  mappertest.hh basisfunctiontest.cc \
mappertest.cc
EXTRA_DIST = 2dgrid.dgf 3dgrid.dgf cube.dgf onesimplex.dgf 

include $(top_srcdir)/am/global-rules

