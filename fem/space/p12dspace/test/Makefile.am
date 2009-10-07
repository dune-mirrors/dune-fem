testdir=$(includedir)/dune/fem/space/p12dspace/test

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

DIMENSION = 2

EXTRADIST = laplace.hh  massoperator.hh  problem.hh  test_l2projection.hh  test_poisson.hh  2dgrid.dgf

if HAVE_DUNE_LOCALFUNCTIONS

check_PROGRAMS          = test_p12dspace 
test_p12dspace_SOURCES  = main.cc ../../../misc/test.cc ../../../misc/suite.cc
test_p12dspace_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DDIMENSION=$(DIMENSION)

endif

include $(top_srcdir)/am/global-rules
