spacedir = $(includedir)/dune/fem/space
space_HEADERS = anisotropicdgspace.hh \
		combineddiscretefunctionspace.hh \
                combinedspace.hh \
                dgspace.hh \
                fvspace.hh \
		genericspace.hh \
                lagrangespace.hh \
                padaptivespace.hh \
                reducedbasisspace.hh

SUBDIRS = anisotropicdgspace basefunctions basisfunctionset \
	  combinedspace common dgspace dofmapper fvspace \
	  genericspace lagrange lagrangespace mapper padaptivespace \
	  reducedbasisspace shapefunctionset test

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = dgspace/libdgspace.la \
                     lagrangespace/liblagrangespace.la

include $(top_srcdir)/am/global-rules
