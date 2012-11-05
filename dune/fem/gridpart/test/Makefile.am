# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=SGRID
GRIDDIM=2

# test programs 
check_PROGRAMS = test-filteredgridpart lagrangeinterpolation test-idgridpart test-geogridpart

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRIDPART -DPOLORDER=1 
test_filteredgridpart_SOURCES = test-filteredgridpart.cc
test_filteredgridpart_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRIDPART

test_idgridpart_SOURCES = test-idgridpart.cc
test_geogridpart_SOURCES = test-geogridpart.cc

EXTRA_DIST = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf \
	     checkgeometry.hh  checkindexset.hh \
	     checkintersections.hh  checkseed.hh \
	     failure.hh

include $(top_srcdir)/am/global-rules
