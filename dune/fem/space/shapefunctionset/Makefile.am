shapefunctionsetdir = $(includedir)/dune/fem/space/shapefunctionset
shapefunctionset_HEADERS = \
	caching.hh \
	lagrange.hh \
	legendre.hh \
	localfunctions.hh \
	mapped.hh \
	orthonormal.hh \
	proxy.hh \
	selectcaching.hh \
	shapefunctionset.hh \
	simple.hh \
	tensorproduct.hh \
	vectorial.hh \
	wrapper.hh

SUBDIRS = orthonormal

EXTRA_DIST = legendrepolynomials.cc

noinst_LTLIBRARIES = libshapefunctionset.la

libshapefunctionset_la_CXXFLAGS = @AM_CPPFLAGS@
libshapefunctionset_la_LIBADD = \
	orthonormal/liborthonormal.la
libshapefunctionset_la_SOURCES = \
	legendrepolynomials.cc
libshapefunctionset_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

include $(top_srcdir)/am/global-rules
