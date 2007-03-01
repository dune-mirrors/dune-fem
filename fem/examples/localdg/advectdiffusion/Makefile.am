LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lcblas
noinst_PROGRAMS = advectdiff dgtest #euler dgtest dgtest1 #datadisp 
advectdiff_SOURCES = main.cc
#euler_SOURCES = euler.cc
dgtest_SOURCES = dgtest.cc
#dgtest1_SOURCES = dgtest.cc
#datadisp_SOURCES = datadisp.cc

GRIDTYPE = YASPGRID
GRIDDIM = 2
PROBLEM = 1
POLORDER = 0

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE) -I/usr/people/robertk/work/Dune/release/dune-istl
CXXFLAGS = -O3 -DNDEBUG -funroll-loops -finline-functions --param max-inline-insns-single=4000 --param large-function-growth=4500 --param inline-unit-growth=4000 -ffast-math -fomit-frame-pointer -march=opteron -msse3 -mfpmath=sse # -Wall
# CXXFLAGS = -pg -O3 -DNDEBUG -funroll-loops -finline-functions --param max-inline-insns-single=4000 --param large-function-growth=4500 --param inline-unit-growth=4000 -ffast-math -march=opteron -msse3 -mfpmath=sse # -Wall
# CXXFLAGS = -g
#-O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500 

