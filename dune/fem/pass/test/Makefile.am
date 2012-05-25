# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

SUBDIRS = dgellipt

# test programs 
check_PROGRAMS = tester 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

tester_SOURCES = main.cc test.cc suite.cc pass_test.cc helper_test.cc

EXTRA_DIST = dgstub.hh \
             helper_test.hh \
             lagrange_fixture.hh \
             passstub.hh \
             pass_test.hh \
             suite.hh \
             test.hh

include $(top_srcdir)/am/global-rules
