POLORDER = 2
#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDDIM = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) ../../libfem.la
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)

ADAPTPROG = lagrangeadapt adapt 
check_PROGRAMS = l2projection lagrangeinterpolation legendretest

l2projection_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
legendretest_SOURCES = legendretest.cc
