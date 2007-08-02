LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lblas

noinst_PROGRAMS = poisson
poisson_SOURCES = poisson.cc

# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
#Default values:
GRIDTYPE=YASPGRID
GRIDDIM=2
POLORDER=1

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -g -Wall --param max-inline-insns-single=1000 --param large-function-growth=1500
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -Wall -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000 --param large-function-growth=1500

