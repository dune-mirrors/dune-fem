POLORDER = 2
WANT_GRAPE = 0
WANT_CACHED = 1

# set default grid type for make all
GRIDTYPE=YASPGRID
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_CONFORM
#GRIDTYPE=ALUGRID_CUBE
GRIDDIM=2

#CXXFLAGS=-g -Wall -Wfatal-errors

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWANT_GRAPE=$(WANT_GRAPE) \
  -DWANT_CACHED_COMM_MANAGER=$(WANT_CACHED)

ADAPTPROG = lagrangeadapt adapt
check_PROGRAMS = $(ADAPTPROG) l2projection lagrangeinterpolation legendretest \
                 lagrangeglobalrefine arraytest dgcomm padapt

l2projection_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
padapt_SOURCES = padapt.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
lagrangeglobalrefine_SOURCES = lagrangeglobalrefine.cc
legendretest_SOURCES = legendretest.cc
arraytest_SOURCES = arraytest.cc
dgcomm_SOURCES = dgcomm.cc

include $(top_srcdir)/am/global-rules

