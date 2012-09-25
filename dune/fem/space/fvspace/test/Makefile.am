GRIDTYPE=YASPGRID
GRIDDIM=2

DIMRANGE = 1

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMRANGE=$(DIMRANGE)

check_PROGRAMS = test-fvspace

test_fvspace_SOURCES = test-fvspace.cc

include $(top_srcdir)/am/global-rules
