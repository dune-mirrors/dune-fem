spacedir = $(includedir)/dune/fem/space

space_HEADERS = anisotropicdgspace.hh \
		combineddiscretefunctionspace.hh \
		combinedspace.hh \
		dgspace.hh \
		discontinuousgalerkin.hh \
		finitevolume.hh \
		fvspace.hh \
		genericspace.hh \
		lagrange.hh \
		lagrangespace.hh \
		padaptivespace.hh \
		reducedbasisspace.hh

SUBDIRS = anisotropicdgspace basefunctions \
	  basisfunctionset combinedspace common \
	  dgspace discontinuousgalerkin dofmapper \
	  finitevolume fvspace genericspace lagrange \
	  lagrangespace mapper padaptivespace \
	  reducedbasisspace shapefunctionset test

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = dgspace/libdgspace.la \
		     discontinuousgalerkin/libdiscontinuousgalerkin.la \
		     lagrange/liblagrange.la \
                     lagrangespace/liblagrangespace.la

include $(top_srcdir)/am/global-rules
