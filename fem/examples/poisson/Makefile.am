LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lblas

noinst_PROGRAMS = poisson
poisson_SOURCES = poisson.cc

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -w -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500

