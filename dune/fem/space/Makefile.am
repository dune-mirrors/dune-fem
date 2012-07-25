spacedir = $(includedir)/dune/fem/space
space_HEADERS = adaptoperator.hh \
                combineddiscretefunctionspace.hh \
                combinedspace.hh \
                dgspace.hh \
                fvspace.hh \
                lagrangespace.hh \
                padaptivespace.hh \
                reducedbasisspace.hh

SUBDIRS = basefunctions common dofmapper mapper test dgspace fvspace \
          lagrangespace combinedspace reducedbasisspace \
          padaptivespace

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = dgspace/libdgspace.la \
                     lagrangespace/liblagrangespace.la

include $(top_srcdir)/am/global-rules
