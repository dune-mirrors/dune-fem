orthonormaldir =  $(includedir)/dune/fem/space/shapefunctionset/orthonormal
orthonormal_HEADERS = \
	orthonormalbase_1d.hh \
	orthonormalbase_2d.hh \
	orthonormalbase_3d.hh

EXTRA_DIST = \
	orthonormalbase_1d.cc \
	orthonormalbase_2d.cc \
	orthonormalbase_3d.cc

noinst_LTLIBRARIES = liborthonormal.la

liborthonormal_la_CXXFLAGS = @AM_CPPFLAGS@
liborthonormal_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB
liborthonormal_la_LIBADD  =
liborthonormal_la_SOURCES = \
	orthonormalbase_1d.cc \
	orthonormalbase_2d.cc \
	orthonormalbase_3d.cc

include $(top_srcdir)/am/global-rules
