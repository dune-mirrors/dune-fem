# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

# test programs 
check_PROGRAMS = test-filteredgridpart lagrangeinterpolation test-idgridpart

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

test_lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
test_lagrangeinterpolation_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRIDPART -DPOLORDER=1 
test_filteredgridpart_SOURCES = test-filteredgridpart.cc
test_filteredgridpart_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRIDPART

test_idgridpart_SOURCES = test-idgridpart.cc

EXTRA_DIST = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf

include $(top_srcdir)/am/global-rules
