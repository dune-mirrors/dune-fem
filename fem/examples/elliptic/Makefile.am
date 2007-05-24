# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
#Default values:
GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER = 1

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lblas

noinst_PROGRAMS = elliptic
elliptic_SOURCES = elliptic.cc

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -DPOLORDER=$(POLORDER) -g --param max-inline-insns-single=1000 --param large-function-growth=1500 -Wall

#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -w -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500 

