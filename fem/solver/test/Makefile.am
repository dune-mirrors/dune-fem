POLORDER = 2
GRIDTYPE = YASPGRID
GRIDDIM = 2

#CXXFLAGS=-O3 -DNOFEL=$(size) -DNDEBUG # -Wall

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)

#check_PROGRAMS = odetest odesystemtest
odetest_SOURCES = odetest.cc
odesystemtest_SOURCES = odesystemtest.cc

include $(top_srcdir)/am/global-rules
