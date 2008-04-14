# IGNORE_HEADERS

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) 
noinst_PROGRAMS = dgtest datadisp 
dgtest_SOURCES = dgtest.cc
datadisp_SOURCES = datadisp.cc

GRIDTYPE = YASPGRID
GRIDDIM = 2
PROBLEM = 1
POLORDER = 1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DGRIDDIM=$(GRIDDIM) -DPROBLEM=$(PROBLEM) -DPOLORDER=$(POLORDER) -D$(GRIDTYPE)

# CXXFLAGS = -O3 -funroll-loops -finline-functions --param max-inline-insns-single=4000 --param large-function-growth=4500 --param inline-unit-growth=4000 -ffast-math -fomit-frame-pointer -march=opteron -msse3 -mfpmath=sse # -Wall
CXXFLAGS = -g -Wall
