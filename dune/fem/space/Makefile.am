spacedir = $(includedir)/dune/fem/space

space_HEADERS = anisotropicdgspace.hh \
		combineddiscretefunctionspace.hh \
		combinedspace.hh \
		discontinuousgalerkin.hh \
		finitevolume.hh \
		genericspace.hh \
		lagrange.hh \
		padaptivespace.hh \
		reducedbasisspace.hh

SUBDIRS = anisotropicdgspace basefunctions \
	  basisfunctionset combinedspace common \
	  discontinuousgalerkin dofmapper \
	  finitevolume genericspace lagrange \
	  mapper padaptivespace reducedbasisspace \
	  shapefunctionset test

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = discontinuousgalerkin/libdiscontinuousgalerkin.la \
		     lagrange/liblagrange.la

include $(top_srcdir)/am/global-rules
