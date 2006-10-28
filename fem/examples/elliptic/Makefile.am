LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lblas

noinst_PROGRAMS = elliptic
elliptic_SOURCES = elliptic.cc

# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  
#Default values:
GRIDTYPE=YASPGRID
GRIDDIM=2

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -w -g --param max-inline-insns-single=1000  --param large-function-growth=1500 

#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM) -w -O3 -unroll -DNDEBUG --param max-inline-insns-single=1000  --param large-function-growth=1500 

