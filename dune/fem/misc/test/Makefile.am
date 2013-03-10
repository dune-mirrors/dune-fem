LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

SUBDIRS =

check_PROGRAMS = example longexample

example_SOURCES = example.cc
longexample_SOURCES = longexample.cc

EXTRA_DIST =
CLEANFILES =

include $(top_srcdir)/am/global-rules
