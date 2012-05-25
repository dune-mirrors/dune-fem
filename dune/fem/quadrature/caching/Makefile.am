SUBDIRS = test

cachingdir =  $(includedir)/dune/fem/quadrature/caching
caching_HEADERS = cacheprovider.hh  pointmapper.hh  pointprovider.hh twistprovider.hh \
cacheprovider.cc  pointprovider.cc  twistprovider.cc twistutility.hh

noinst_LTLIBRARIES = libcaching.la

libcaching_la_SOURCES = twistutility.cc
libcaching_la_CPPFLAGS = $(UG_CPPFLAGS) $(AM_CPPFLAGS) -DCOMPILE_FEMLIB
libcaching_la_LIBADD = # $(UG_LDFLAGS) $(UG_LIBS)

include $(top_srcdir)/am/global-rules
