LDADD=$(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_LDFLAGS= $(ALL_PKG_LDFLAGS)
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS)

TESTS = odetest test-timeprovider
check_PROGRAMS = $(TESTS)

odetest_SOURCES = odetest.cc
test_timeprovider_SOURCES = test-timeprovider.cc

EXTRA_DIST = parameter

include $(top_srcdir)/am/global-rules
