commtestdir = $(includedir)/dune/fem/petsc/experiments/commtest

#GRIDTYPE = ALUGRID_SIMPLEX
#GRIDTYPE = ALUGRID_CUBE
#GRIDTYPE = SPGRID
GRIDTYPE = YASPGRID
#YASPGRID = 1
GRIDDIM = 2
DIMRANGE = 2
POLORDER = 1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DFEMTIMER -DGRIDDIM=$(GRIDDIM) -DDIMRANGE=$(DIMRANGE) \
			  -DPOLORDER=$(POLORDER) 
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)

CHECK_PROGS = 

check_PROGRAMS = $(CHECK_PROGS)

CLEANFILES = *.tex data/*

# commtest_SOURCES = commtest.cc
# commtest_CPPFLAGS = $(AM_CPPFLAGS)

EXTRA_DIST = domain_commtest.dgf parameter_commtest

include $(top_srcdir)/am/global-rules
