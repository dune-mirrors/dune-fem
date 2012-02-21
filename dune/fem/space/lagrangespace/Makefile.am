AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = test 

lagrangespacedir = $(includedir)/dune/fem/space/lagrangespace
lagrangespace_HEADERS = adaptmanager.hh  basefunctions.hh \
                        basefunctions_inline.hh  dofmappercode.hh \
                        genericbasefunctions.hh \
                        genericgeometry.hh  genericlagrangepoints.hh \
                        lagrangepoints.hh \
                        lagrangespace.hh  mapper.hh  restrictprolong.hh

noinst_LTLIBRARIES = liblagrangespace.la

liblagrangespace_la_SOURCES = basefunctions.cc
liblagrangespace_la_CPPFLAGS = $(AM_CPPFLAGS)
liblagrangespace_la_LIBADD =

include $(top_srcdir)/am/global-rules

