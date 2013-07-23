POLORDER = 2
DIMRANGE = 3
WANT_GRAPE = 0

# set default grid type for make all
GRIDTYPE=YASPGRID
#GRIDTYPE=PARALLELGRID_ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_SIMPLEX
#GRIDTYPE=ALUGRID_CONFORM
#GRIDTYPE=ALUGRID_CUBE
#GRIDTYPE=ALBERTAGRID
GRIDDIM=3

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DDIMRANGE=$(DIMRANGE) -DWANT_GRAPE=$(WANT_GRAPE)

ADAPTPROG = lagrangeadapt adapt padapt

# programs run during make check 
TESTS = $(ADAPTPROG) dgcomm-4.sh dgcomm-cached-4.sh

check_PROGRAMS = $(ADAPTPROG) \
		 arraytest \
		 combinedspacetest \
		 dgcomm \
		 dgcomm_cached \
		 interpolation \
     l2projection \
     lagrangeglobalrefine \
		 lagrangeinterpolation \
		 legendretest

EXTRA_DIST = systemmatrix.cc dgcomm-4.sh dgcomm-cached-4.sh parameter *.dgf

l2projection_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
padapt_SOURCES = padapt.cc
combinedspacetest_SOURCES = combinedspacetest.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
interpolation_SOURCES = interpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
lagrangeglobalrefine_SOURCES = lagrangeglobalrefine.cc
legendretest_SOURCES = legendretest.cc
arraytest_SOURCES = arraytest.cc
dgcomm_SOURCES = dgcomm.cc
dgcomm_cached_SOURCES = dgcomm.cc
dgcomm_CPPFLAGS = $(AM_CPPFLAGS) -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3 -DWANT_CACHED_COMM_MANAGER=0
dgcomm_cached_CPPFLAGS = $(AM_CPPFLAGS) -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3 

include $(top_srcdir)/am/global-rules
