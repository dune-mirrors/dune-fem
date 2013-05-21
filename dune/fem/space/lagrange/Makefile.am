AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = test 

lagrangedir = $(includedir)/dune/fem/space/lagrange

lagrange_HEADERS = adaptmanager.hh \
		   dofmappercode.hh \
		   genericlagrangepoints.hh \
		   shapefunctionset.hh \
		   capabilities.hh \
		   genericbasefunctions.hh \
		   lagrangepoints.hh \
		   space.hh \
		   declaration.hh \
		   genericgeometry.hh \
		   restrictprolong.hh \
		   storage.hh

noinst_LTLIBRARIES = liblagrange.la

liblagrange_la_SOURCES = shapefunctionset.cc
liblagrange_la_CPPFLAGS = $(AM_CPPFLAGS)
liblagrange_la_LIBADD =

include $(top_srcdir)/am/global-rules
