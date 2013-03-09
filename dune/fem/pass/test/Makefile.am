GRIDTYPE = YASPGRID
GRIDDIM = 2

DIMRANGE = 2
POLORDER = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

SUBDIRS =

# test programs 
check_PROGRAMS = test-pass

test_pass_SOURCES = test-pass.cc
test_pass_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE) -DPOLORDER=$(POLORDER)

DGFFILES = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf

EXTRA_DIST = $(DGFFILES) parameter

CLEANFILES =

include $(top_srcdir)/am/global-rules
