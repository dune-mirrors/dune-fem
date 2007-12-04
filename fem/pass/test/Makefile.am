SUBDIRS = dgellipt

# test programs 
check_PROGRAMS = tester 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

# check_PROGRAMS = tester 
tester_SOURCES = main.cc test.cc suite.cc pass_test.cc helper_test.cc \
	test.hh suite.hh pass_test.hh helper_test.hh passstub.hh dgstub.hh
