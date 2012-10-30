AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = test 

lagrangedir = $(includedir)/dune/fem/space/lagrange

lagrange_HEADERS = adaptmanager.hh  declaration.hh \
		   dofmappercode.hh  genericbasefunctions.hh \
		   genericgeometry.hh  genericlagrangepoints.hh \
		   lagrangepoints.hh  lagrangespace.hh \
		   restrictprolong.hh  shapefunctionset.hh \
		   space.hh  storage.hh

noinst_LTLIBRARIES = liblagrange.la

liblagrange_la_SOURCES = shapefunctionset.cc
liblagrange_la_CPPFLAGS = $(AM_CPPFLAGS)
liblagrange_la_LIBADD =

include $(top_srcdir)/am/global-rules
