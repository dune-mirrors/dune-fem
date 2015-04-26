POLORDER = 2
DIMRANGE = 3

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DDIMRANGE=$(DIMRANGE)

DFPROG = dftest

# programs run during make check 
TESTS = $(DFPROG)

check_PROGRAMS = $(DFPROG)
dftest_SOURCES = dftest.cc

include $(top_srcdir)/am/global-rules
