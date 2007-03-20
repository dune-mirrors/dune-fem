# $Id$

GRIDTYPE=YASPGRID
#GRIDTYPE=ALBERTAGRID
#GRIDTYPE=ALUGRID_SIMPLEX
CXXFLAGS=-O3 -funroll-loops -DNDEBUG -Wall
DGELLIPTPROG = dgellipt 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS) -I$(HOME)/work/Dune/release/dune-istl $(MPI_CPPFLAGS)

check_PROGRAMS = $(DGELLIPTPROG) 
dgellipt_SOURCES = main.cc

EXTRA_DIST = problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log
