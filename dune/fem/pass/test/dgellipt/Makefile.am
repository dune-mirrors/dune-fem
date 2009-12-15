# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

#GRIDTYPE=ALBERTAGRID
#GRIDTYPE=ALUGRID_SIMPLEX
#WORLDDIM=2
#GRIDTYPE=ALUGRID_CUBE
#CXXFLAGS=-O3 -funroll-loops -DNDEBUG -Wall
#CXXFLAGS = -g -w
#CXXFLAGS = -DGRIDTYPE=ALBERTAGRID -g -O2 -Wall
DGELLIPTPROG = dgellipt 
#CXX = mpiCC

#LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS) $(MPI_CPPFLAGS)

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNE_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = $(DGELLIPTPROG) 
dgellipt_SOURCES = main.cc

EXTRA_DIST = benchmark.cc problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log

# CXXFLAGS =

include $(top_srcdir)/am/global-rules

