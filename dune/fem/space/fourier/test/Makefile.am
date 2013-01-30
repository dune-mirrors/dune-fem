GRIDTYPE=YASPGRID
GRIDDIM=2

ORDER = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DORDER=$(ORDER)

check_PROGRAMS = test-space

test_space_SOURCES = test-space.cc

EXTRA_DIST = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf

CLEANFILES = *.vtu *.pvtu *.log

include $(top_srcdir)/am/global-rules
