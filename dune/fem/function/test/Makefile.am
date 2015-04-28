# set default grid type for make all
GRIDTYPE = YASPGRID
GRIDDIM = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

DFPROG = dftest

# programs run during make check 
TESTS = $(DFPROG)

check_PROGRAMS = $(DFPROG)
dftest_SOURCES = dftest.cc

EXTRA_DIST = \
	1dgrid.dgf \
	2dgrid.dgf \
	3dgrid.dgf

include $(top_srcdir)/am/global-rules
