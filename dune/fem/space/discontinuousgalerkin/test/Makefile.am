GRIDTYPE=SPGRID
GRIDDIM=2

POLORDER = 2
DIMRANGE = 1

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DDIMRANGE=$(DIMRANGE)

check_PROGRAMS = test-dgspaces

test_dgspaces_SOURCES = test-dgspaces.cc

CLEANFILES = *.gnu *.vtu *.pvtu *.vtp *.pvtp *.log 

include $(top_srcdir)/am/global-rules
