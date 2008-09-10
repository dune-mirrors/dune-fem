AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = test 

lagrangespacedir = $(includedir)/dune/fem/space/lagrangespace
lagrangespace_HEADERS = adaptmanager.hh \
			basefunctions.hh \
			basefunctions_inline.hh \
			genericbasefunctions.hh \
			genericgeometry.hh \
			genericlagrangepoints.hh \
			lagrangedatahandle.hh \
			lagrangepoints.hh \
			lagrangespace.hh \
			mapper.hh

noinst_LTLIBRARIES = liblagrangespace.la

liblagrangespace_la_SOURCES = basefunctions.cc
liblagrangespace_la_CPPFLAGS = $(AM_CPPFLAGS) -DCOMPILE_FEMLIB
liblagrangespace_la_LIBADD =

include $(top_srcdir)/am/global-rules

