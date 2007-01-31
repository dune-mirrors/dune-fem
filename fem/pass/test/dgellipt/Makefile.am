# $Id$

GRIDTYPE=ALBERTAGRID
GRIDDIM=2

if ALBERTA
DGELLIPTPROG = dgellipt 
endif 
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS)

noinst_PROGRAMS = $(DGELLIPTPROG) 
if ALBERTA
dgellipt_SOURCES = main.cc
endif

EXTRA_DIST = problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log
