# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

# test programs 
check_PROGRAMS = test-filteredgridpart 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

test_filteredgridpart_SOURCES = test-filteredgridpart.cc
test_filteredgridpart_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRIDPART

EXTRA_DIST = radialfilter.hh 

include $(top_srcdir)/am/global-rules
