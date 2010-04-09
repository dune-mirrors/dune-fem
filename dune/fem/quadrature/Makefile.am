SUBDIRS = caching test

quadraturedir =  $(includedir)/dune/fem/quadrature
quadrature_HEADERS = cachingpointlist.hh \
                     cachingquadrature.hh \
                     defaultquadratures.hh \
                     dunequadratures.hh \
                     elementpointlist.hh \
                     elementpointlistbase.hh \
                     elementquadrature.hh \
                     femquadratures.hh \
                     femquadratures_inline.hh \
                     gausspoints.hh \
                     gausspoints_implementation.hh \
                     idprovider.hh \
                     integrator.hh \
                     quadrature.hh \
                     quadratureimp.hh \
                     quadratureimp_inline.hh \
                     pyramidpoints.hh \
                     quadprovider.hh \
                     simplexpoints.hh \
                     intersectionquadrature.hh

noinst_LTLIBRARIES = libquadrature.la

libquadrature_la_SOURCES =
libquadrature_la_LIBADD = caching/libcaching.la

include $(top_srcdir)/am/global-rules
