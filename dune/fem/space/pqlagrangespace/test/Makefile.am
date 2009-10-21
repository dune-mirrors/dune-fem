# IGNORE_HEADERS

testdir=$(includedir)/dune/fem/space/pqlagrangespace/test

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) -lgmpxx -lgmp -lmpfr
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

EXTRA_DIST = laplace.hh  massoperator.hh  problem.hh  \
             test_l2projection.hh  test_poisson.hh \
             2dgrid.dgf 3dgrid.dgf
#test_HEADERS = laplace.hh test_l2projection.hh \
#               massoperator.hh  test_poisson.hh \
#               problem.hh


if HAVE_DUNE_LOCALFUNCTIONS

check_PROGRAMS          = test_pqlagrangespace 
test_pqlagrangespace_SOURCES  = main.cc ../../../misc/test.cc ../../../misc/suite.cc
test_pqlagrangespace_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DHAVE_GMP=0 -DHAVE_ALGLIB=1

endif

include $(top_srcdir)/am/global-rules
