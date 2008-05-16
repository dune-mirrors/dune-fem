# IGNORE_HEADERS

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) -lblas

POISSON_CPPFLAGS = -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -DPROBLEM=$(PROBLEM)
OPTIM_CPPFLAGS = -g -Wall --param max-inline-insns-single=1000 --param large-function-growth=1500



noinst_PROGRAMS = serial parallel
serial_SOURCES = poisson.cc
serial_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(OPTIM_CPPFLAGS) -DUSE_GRAPE=1

parallel_SOURCES = poisson.cc
parallel_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(OPTIM_CPPFLAGS) $(MPI_CPPFLAGS)
# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
#Default values:
#GRIDTYPE=SGRID
GRIDTYPE=YASPGRID
GRIDDIM=2
POLORDER=1
PROBLEM=1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(POISSON_CPPFLAGS) $(OPTIM_CPPFLAGS)
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -Wall -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000 --param large-function-growth=1500

