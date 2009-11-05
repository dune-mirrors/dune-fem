# IGNORE_HEADERS

POLORDER=3

testdir=$(includedir)/dune/fem/space/pqlagrangespace/test

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) 
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

EXTRA_DIST = laplace.hh  massoperator.hh  problem.hh  \
             test_l2projection.hh  test_poisson.hh \
             2dgrid.dgf 3dgrid.dgf hybrid2_0.dgf
#test_HEADERS = laplace.hh test_l2projection.hh \
#               massoperator.hh  test_poisson.hh \
#               problem.hh


if HAVE_DUNE_LOCALFUNCTIONS

check_PROGRAMS          = test_pqlagrangespace test-pk2d

test_pqlagrangespace_SOURCES  = main.cc ../../../misc/test.cc ../../../misc/suite.cc
test_pqlagrangespace_CPPFLAGS = $(ALL_PKG_CPPFLAGS) 

test_pk2d_SOURCES = test-pk2d.cc
test_pk2d_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER)

endif

include $(top_srcdir)/am/global-rules
