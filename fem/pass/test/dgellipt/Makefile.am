# $Id$

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS) $(MPI_CPPFLAGS)

noinst_PROGRAMS = dgellipt 
dgellipt_SOURCES = main.cc

CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log
