# $Id$

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS= -w -O3 -unroll -DNDEBUG

noinst_PROGRAMS = advectdiff euler advection #datadisp
advectdiff_SOURCES = main.cc
euler_SOURCES = euler.cc
advection_SOURCES = advection.cc
#datadisp_SOURCES = datadisp.cc

