POLORD=0

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -w -DPOLORDER=$(POLORD)

noinst_PROGRAMS = l2pro adapt
l2pro_SOURCES = l2projection.cc
adapt_SOURCES = adapt.cc
