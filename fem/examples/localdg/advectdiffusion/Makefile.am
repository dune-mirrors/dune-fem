LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

noinst_PROGRAMS = advectdiff dgtest #euler dgtest dgtest1 #datadisp 
advectdiff_SOURCES = main.cc
#euler_SOURCES = euler.cc
dgtest_SOURCES = dgtest.cc
#dgtest1_SOURCES = dgtest.cc
#datadisp_SOURCES = datadisp.cc

GRIDTYPE = YASPGRID
DUNE_WORLD_DIM = $(GRIDDIM)
PROBLEM = 1
POLORDER = 0

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE) -DDUNE_PROBLEM_DIM=$(DUNE_PROBLEM_DIM) -w -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500 -DGRIDDIM=$(DUNE_PROBLEM_DIM) -DGRIDDIM=$(GRIDDIM)

