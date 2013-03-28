GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER = 2
DIMRANGE = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = test-vectorialbasisfunctionset

test_vectorialbasisfunctionset_SOURCES = test-vectorialbasisfunctionset.cc
test_vectorialbasisfunctionset_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE) -DPOLORDER=$(POLORDER)

DGFFILES = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf

EXTRA_DIST = $(DGFFILES)

CLEANFILES = dgfparser.log

include $(top_srcdir)/am/global-rules
