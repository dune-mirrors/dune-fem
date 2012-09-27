DIMRANGE = 3
POLORDER = 2

GRIDTYPE=YASPGRID
GRIDDIM=3

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMRANGE=$(DIMRANGE) -DPOLORDER=$(POLORDER)

check_PROGRAMS = test-multiindexset test-space

test_multiindexset_SOURCES = test-multiindexset.cc
test_space_SOURCES = test-space.cc

EXTRA_DIST =

include $(top_srcdir)/am/global-rules
