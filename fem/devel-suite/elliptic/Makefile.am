# do not change these values, but use 
# make GRIDTYPE=ALUGRID_SIMPLEX GRIDDIM=3
# etc. for overriding these settings  

# IGNORE_HEADERS

#Default values:
GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER = 1
PROBLEM = POISSON

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)

noinst_PROGRAMS = elliptic
elliptic_SOURCES = elliptic.cc

AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -D$(PROBLEM) -Wall -DNDEBUG -O3 --param max-inline-insns-single=1000 --param large-function-growth=1500
#AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -Wall -g --param max-inline-insns-single=1000 --param large-function-growth=1500
