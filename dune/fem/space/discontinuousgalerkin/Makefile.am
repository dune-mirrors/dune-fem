noinst_LTLIBRARIES = libdiscontinuousgalerkin.la

libdiscontinuousgalerkin_la_CXXFLAGS = @AM_CPPFLAGS@
libdiscontinuousgalerkin_la_LIBADD  =

libdiscontinuousgalerkin_la_SOURCES = dgbasefunctions_1.cc \
                        dgbasefunctions_2.cc \
                        dgbasefunctions_3.cc \
			legendrepoly.cc
libdiscontinuousgalerkin_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

discontinuousgalerkindir =  $(includedir)/dune/fem/space/discontinuousgalerkin
discontinuousgalerkin_HEADERS = declaration.hh  default.hh  dgbasefunctions.hh \
				dgspace.hh  lagrange.hh  legendredgbasefunctions.hh \
				legendre.hh  legendrepoly.hh  localdgmassmatrix.hh \
				localrestrictprolong.hh  orthonormalbase_mod.hh  space.hh

EXTRA_DIST = orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc

include $(top_srcdir)/am/global-rules

