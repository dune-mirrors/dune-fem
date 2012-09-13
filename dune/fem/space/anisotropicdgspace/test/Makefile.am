POLORDER = 2

GRIDTYPE=YASPGRID
GRIDDIM=3

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)

check_PROGRAMS = test-multiindexset

test_multiindexset_SOURCES = test-multiindexset.cc

EXTRA_DIST =

include $(top_srcdir)/am/global-rules
