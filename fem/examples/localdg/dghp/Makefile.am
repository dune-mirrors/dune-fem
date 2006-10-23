LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CXXFLAGS= 
CXXFLAGS=-O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500 
# CXXFLAGS=-g
AM_LDFLAGS=

noinst_PROGRAMS = dgtest datadisp mhd3dvis
dgtest_SOURCES = dgtest.cc
datadisp_SOURCES = datadisp.cc
mhd3dvis_SOURCES = mhd3dvis.cc

GRIDTYPE = DEFAULTGRID
DUNE_WORLD_DIM = $(DUNE_PROBLEM_DIM)
PROBLEM = 1
POLORDER = 0

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DGRIDDIMWORLD=$(GRIDDIM)

