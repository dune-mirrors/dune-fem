noinst_LTLIBRARIES = libdiscontinuousgalerkin.la

libdiscontinuousgalerkin_la_CXXFLAGS = @AM_CPPFLAGS@
libdiscontinuousgalerkin_la_LIBADD  =

libdiscontinuousgalerkin_la_SOURCES = dgbasefunctions_1.cc \
                        dgbasefunctions_2.cc \
                        dgbasefunctions_3.cc \
			legendrepolynomials.cc
libdiscontinuousgalerkin_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

discontinuousgalerkindir =  $(includedir)/dune/fem/space/discontinuousgalerkin
discontinuousgalerkin_HEADERS = declaration.hh \
				lagrange.hh \
				localdgmassmatrix.hh \
				shapefunctionset.hh \
				default.hh \
				legendre.hh \
				localrestrictprolong.hh \
				space.hh \
				dgbasefunctions.hh \
				legendrepolynomials.hh \
				orthonormalbase_mod.hh \

EXTRA_DIST = orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc

include $(top_srcdir)/am/global-rules
