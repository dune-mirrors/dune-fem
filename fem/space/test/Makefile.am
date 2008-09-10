POLORDER = 2
WANT_GRAPE = 0
#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDDIM = 2

CXXFLAGS=-g -Wall

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWANT_GRAPE=$(WANT_GRAPE)

ADAPTPROG = lagrangeadapt adapt 
check_PROGRAMS = $(ADAPTPROG) l2projection lagrangeinterpolation legendretest combinedtest

l2projection_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
combinedtest_SOURCES = combinedtest.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
legendretest_SOURCES = legendretest.cc

include $(top_srcdir)/am/global-rules

