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

ADAPTPROG = adapt adapt_combinedspace_pointbased adapt_combinedspace_variablebased lagrangeadapt functiontupleadapt padapt

# programs run during make check 
TESTS = $(ADAPTPROG) dgcomm-4.sh dgcomm-cached-4.sh

check_PROGRAMS = $(ADAPTPROG) \
		 adapt \
		 adapt_combinedspace_pointbased \
		 adapt_combinedspace_variablebased \
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
adapt_combinedspace_pointbased_SOURCES = adapt.cc
adapt_combinedspace_variablebased_SOURCES = adapt.cc
padapt_SOURCES = padapt.cc
combinedspacetest_SOURCES = combinedspacetest.cc
lagrangeinterpolation_SOURCES = lagrangeinterpolation.cc
interpolation_SOURCES = interpolation.cc
lagrangeadapt_SOURCES = lagrangeadapt.cc
functiontupleadapt_SOURCES = functiontupleadapt.cc
lagrangeglobalrefine_SOURCES = lagrangeglobalrefine.cc
legendretest_SOURCES = legendretest.cc
arraytest_SOURCES = arraytest.cc
dgcomm_SOURCES = dgcomm.cc
dgcomm_cached_SOURCES = dgcomm.cc
dgcomm_CPPFLAGS = $(AM_CPPFLAGS) -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3 -DWANT_CACHED_COMM_MANAGER=0
dgcomm_cached_CPPFLAGS = $(AM_CPPFLAGS) -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3 
adapt_combinedspace_pointbased_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_COMBINED_SPACE -DPOINTBASED
adapt_combinedspace_variablebased_CPPFLAGS = $(AM_CPPFLAGS) -DUSE_COMBINED_SPACE

include $(top_srcdir)/am/global-rules
