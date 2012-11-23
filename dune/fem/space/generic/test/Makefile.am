GRIDTYPE=YASPGRID
GRIDDIM=2

DIMRANGE = 3

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMRANGE=$(DIMRANGE)

check_PROGRAMS = test-space

test_space_SOURCES = test-space.cc

EXTRA_DIST = exactsolution.hh

CLEANFILES = *.vtu *.pvtu *.log

include $(top_srcdir)/am/global-rules
