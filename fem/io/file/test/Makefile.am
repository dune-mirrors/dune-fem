POLORDER = 2
WANT_GRAPE = 0

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

check_PROGRAMS = persistencemanagertest # dataoutputtest

persistencemanagertest_SOURCES = persistencemanagertest.cc 
dataoutputtest_SOURCES = dataoutputtest.cc
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -DPOLORDER=$(POLORDER) -DWANT_GRAPE=$(WANT_GRAPE)

CXXFLAGS = -g

include $(top_srcdir)/am/global-rules

