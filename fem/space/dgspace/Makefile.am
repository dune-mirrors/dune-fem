noinst_LTLIBRARIES = libdgspace.la

libdgspace_la_CXXFLAGS = @AM_CPPFLAGS@
libdgspace_la_LIBADD  =

libdgspace_la_SOURCES = dgbasefunctions_1.cc \
                        dgbasefunctions_2.cc \
			dgbasefunctions_3.cc
libdgspace_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_LIB

dgspacedir =  $(includedir)/dune/fem/space/dgspace
dgspace_HEADERS = dgbasefunctions.hh  legendredgbasefunctions.hh \
dgadaptiveleafgridpart.hh  dgleafindexset.hh   legendrepoly.hh \
dgadaptiveleafindexset.hh  dgmapper.hh orthonormalbase_mod.hh \
dgadaptmanager.hh dgspace.hh dgadaptoperator.hh dgdatahandle.hh
