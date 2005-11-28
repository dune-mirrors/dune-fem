# $Id$

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -L/hosts/aragorn/raid5/mario/Dune/dune/lib -ldune
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS= -w -O3 -unroll -DNDEBUG

noinst_PROGRAMS = transport euler advection
transport_SOURCES = main.cc
euler_SOURCES = euler.cc
advection_SOURCES = advection.cc

