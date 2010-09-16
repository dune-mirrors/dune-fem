# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

# also distribute the headers for the checks
EXTRA_DIST = cache_test.hh  refelem_test.hh  twist_test.hh  point_test.hh

check_PROGRAMS = tester 
tester_SOURCES = main.cc ../../../misc/test.cc ../../../misc/suite.cc \
	refelem_test.cc cache_test.cc twist_test.cc point_test.cc

include $(top_srcdir)/am/global-rules

