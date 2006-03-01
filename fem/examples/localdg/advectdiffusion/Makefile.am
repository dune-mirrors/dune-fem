LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)

noinst_PROGRAMS = advectdiff euler datadisp dgtest
advectdiff_SOURCES = main.cc
euler_SOURCES = euler.cc
dgtest_SOURCES = dgtest.cc
datadisp_SOURCES = datadisp.cc

GRIDTYPE = DEFAULTGRID
DUNE_WORLD_DIM = $(DUNE_PROBLEM_DIM)
PROBLEM = 1
POLORDER = 0

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE) -DGRIDDIMENSION=$(DUNE_PROBLEM_DIM) -w -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500

