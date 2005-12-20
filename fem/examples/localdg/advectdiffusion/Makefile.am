# $Id$

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -L/hosts/aragorn/raid5/mario/Dune/dune/lib -ldune
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS= -w -O3 -unroll -DNDEBUG

noinst_PROGRAMS = advectdiff euler shallow advection datadisp
advectdiff_SOURCES = main.cc
euler_SOURCES = euler.cc
shallow_SOURCES = shallow.cc
advection_SOURCES = advection.cc
datadisp_SOURCES = datadisp.cc

