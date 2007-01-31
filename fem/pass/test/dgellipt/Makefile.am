# $Id$

#GRIDTYPE=YASPGRID
GRIDTYPE=ALBERTAGRID
GRIDDIM=2
#GRIDTYPE=ALUGRID_CUBE
#GRIDTYPE=ALUGRID_SIMPLEX
#CXXFLAGS=-O3 -funroll-loops -DNDEBUG -pg 
if ALBERTA
DGELLIPTPROG = dgellipt 
endif 
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS) $(MPI_CPPFLAGS)

noinst_PROGRAMS = $(DGELLIPTPROG) 
if ALBERTA
dgellipt_SOURCES = main.cc
endif

EXTRA_DIST = problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log
