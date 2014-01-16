discontinuousgalerkindir =  $(includedir)/dune/fem/space/discontinuousgalerkin
discontinuousgalerkin_HEADERS = capabilities.hh \
				declaration.hh \
				default.hh \
				hierarchiclegendre.hh	\
				hierarchiclegendremap.hh	\
				lagrange.hh \
				legendre.hh \
				legendrepolynomials.hh \
				localdgmassmatrix.hh \
				localinterpolation.hh \
				localrestrictprolong.hh \
				orthonormalbase_1d.hh \
				orthonormalbase_2d.hh \
				orthonormalbase_3d.hh \
				orthonormalbase_mod.hh \
				shapefunctionset.hh \
				space.hh

EXTRA_DIST = orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc

noinst_LTLIBRARIES = libdiscontinuousgalerkin.la

libdiscontinuousgalerkin_la_CXXFLAGS = @AM_CPPFLAGS@
libdiscontinuousgalerkin_la_LIBADD  =
libdiscontinuousgalerkin_la_SOURCES = legendrepolynomials.cc \
				      orthonormalbase_1d.cc \
				      orthonormalbase_2d.cc \
				      orthonormalbase_3d.cc
libdiscontinuousgalerkin_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

include $(top_srcdir)/am/global-rules
