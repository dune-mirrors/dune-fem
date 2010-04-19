POLORDER = 2
WANT_GRAPE = 0

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

CXXFLAGS=-g -Wall

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWANT_GRAPE=$(WANT_GRAPE)

ADAPTPROG = lagrangeadapt adapt 
check_PROGRAMS = $(ADAPTPROG) l2projection lagrangeinterpolation legendretest \
                 lagrangeglobalrefine

l2projection_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
lagrangeglobalrefine_SOURCES = lagrangeglobalrefine.cc
legendretest_SOURCES = legendretest.cc

include $(top_srcdir)/am/global-rules

