# $Id$

DGELLIPTPROG = dgellipt 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS)

check_PROGRAMS = $(DGELLIPTPROG) 
dgellipt_SOURCES = main.cc

EXTRA_DIST = problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log
