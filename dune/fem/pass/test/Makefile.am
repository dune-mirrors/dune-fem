AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(DUNEMPILIBS)

check_PROGRAMS = test-insertoperatorpass
TESTS = test-insertoperatorpass

test_insertoperatorpass_SOURCES = test-insertoperatorpass.cc

SUBDIRS =

EXTRA_DIST = CMakeLists.txt

CLEANFILES =

include $(top_srcdir)/am/global-rules
