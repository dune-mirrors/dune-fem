LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

check_PROGRAMS = tester generic 
tester_SOURCES = main.cc ../../../misc/test.cc ../../../misc/suite.cc 
generic_SOURCES = generic.cc

testdir=$(includedir)/fem/space/lagrangespace/test
test_HEADERS = basefunctiontest.hh  mappertest.hh basefunctiontest.cc \
mappertest.cc
EXTRA_DIST = cube.dgf  onesimplex.dgf
