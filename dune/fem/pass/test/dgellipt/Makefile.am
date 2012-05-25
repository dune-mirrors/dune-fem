# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

DGELLIPTPROG = dgellipt 

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = $(DGELLIPTPROG) 
dgellipt_SOURCES = main.cc

EXTRA_DIST = benchmark.cc problem.cc discretemodels.hh  discretization.hh  model.hh
CLEANFILES = dgfparser.log grid2d.dgf.* manager.*.log

include $(top_srcdir)/am/global-rules
