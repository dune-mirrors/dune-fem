LDADD=$(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_LDFLAGS= $(ALL_PKG_LDFLAGS)
AM_CPPFLAGS += $(ALL_PKG_CPPFLAGS)

TESTS = odetest
check_PROGRAMS = odetest

odetest_SOURCES = odetest.cc

include $(top_srcdir)/am/global-rules
