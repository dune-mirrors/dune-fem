# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = dgspace/libdgspace.la

SUBDIRS = test common test dgspace fvspace lagrangespace combinedspace \
basefunctions reducedbasisspace 

spacedir =  $(includedir)/dune/fem/space
space_HEADERS = dgspace.hh lagrangespace.hh combinedspace.hh fvspace.hh \
reducedbasisspace.hh
