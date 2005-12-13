LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

check_PROGRAMS = tester 
tester_SOURCES = main.cc ../../misc/test.cc ../../misc/suite.cc \
	refelem_test.cc point_test.cc twist_test.cc cache_test.cc \
	../../quadrature/idprovider.cc ../../quadrature/gausspoints.cc 
