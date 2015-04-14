GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER = 3
DIMRANGE = 5

USE_VERTICAL_DOF_ALIGNMENT = 1

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = test-vectorialbasisfunctionset test-defaultbasisfunctionset test-tuplebasisfunctionset test-simplebasisfunctionset
TESTS = test-vectorialbasisfunctionset test-defaultbasisfunctionset test-tuplebasisfunctionset test-simplebasisfunctionset

test_vectorialbasisfunctionset_SOURCES = test-vectorialbasisfunctionset.cc
test_vectorialbasisfunctionset_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE) -DPOLORDER=$(POLORDER) -DUSE_VERTICAL_DOF_ALIGNMENT=$(USE_VERTICAL_DOF_ALIGNMENT)

test_defaultbasisfunctionset_SOURCES = test-defaultbasisfunctionset.cc
test_defaultbasisfunctionset_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE)

test_tuplebasisfunctionset_SOURCES = test-tuplebasisfunctionset.cc
test_tuplebasisfunctionset_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE)

test_simplebasisfunctionset_SOURCES = test-simplebasisfunctionset.cc
test_simplebasisfunctionset_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DDIMRANGE=$(DIMRANGE)

DGFFILES = 1dgrid.dgf 2dgrid.dgf 3dgrid.dgf

EXTRA_DIST = $(DGFFILES)

CLEANFILES = dgfparser.log

basisfunctionsettestdir = $(includedir)/dune/fem/space/basisfunctionset/test
basisfunctionsettest_HEADERS = checkbasisfunctionset.hh

include $(top_srcdir)/am/global-rules
