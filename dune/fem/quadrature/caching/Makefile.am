SUBDIRS = test

cachingdir =  $(includedir)/dune/fem/quadrature/caching
caching_HEADERS = cacheprovider.cc  pointmapper.hh    pointprovider.hh  topology.cc  twistprovider.cc  twistutility.cc \
									cacheprovider.hh  pointprovider.cc  registry.hh       topology.hh  twistprovider.hh  twistutility.hh


noinst_LTLIBRARIES = libcaching.la

libcaching_la_SOURCES = twistutility.cc topology.cc
libcaching_la_CPPFLAGS = $(UG_CPPFLAGS) $(AM_CPPFLAGS) -DCOMPILE_FEMLIB
libcaching_la_LIBADD = # $(UG_LDFLAGS) $(UG_LIBS)

include $(top_srcdir)/am/global-rules
