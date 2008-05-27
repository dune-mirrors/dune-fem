noinst_LTLIBRARIES = libdgspace.la

libdgspace_la_CXXFLAGS = @AM_CPPFLAGS@
libdgspace_la_LIBADD  =

libdgspace_la_SOURCES = dgbasefunctions_1.cc \
                        dgbasefunctions_2.cc \
			dgbasefunctions_3.cc
libdgspace_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

dgspacedir =  $(includedir)/dune/fem/space/dgspace
dgspace_HEADERS = dgadaptiveleafgridpart.hh \
                  dgadaptiveleafindexset.hh \
                  dgadaptmanager.hh \
		  dgadaptoperator.hh \
		  dgbasefunctions.hh \
		  dgdatahandle.hh \
		  dgleafindexset.hh \
		  dgmapper.hh \
		  dgspace.hh \
                  legendredgbasefunctions.hh \
		  legendrepoly.hh \
		  orthonormalbase_mod.hh

EXTRA_DIST = orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc
