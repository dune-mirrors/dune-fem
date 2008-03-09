POLORDER = 2
GRIDTYPE = YASPGRID
GRIDDIM = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) ../../libfem.la
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)

# check_PROGRAMS = odetest odesystemtest

odetest_SOURCES = odetest.cc
odesystemtest_SOURCES = odesystemtest.cc
CXXFLAGS=-O3 -DNOFEL=$(size) -DNDEBUG # -Wall
