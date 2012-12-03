spacedir = $(includedir)/dune/fem/space

space_HEADERS = combineddiscretefunctionspace.hh \
		combinedspace.hh \
		discontinuousgalerkin.hh \
		finitevolume.hh \
		lagrange.hh \
		padaptivespace.hh \
		rannacherturek.hh \
		reducedbasisspace.hh

SUBDIRS = basefunctions basisfunctionset \
	  combinedspace common \
	  discontinuousgalerkin dofmapper \
	  finitevolume lagrange mapper \
	  padaptivespace rannacherturek \
	  reducedbasisspace shapefunctionset \
	  test

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = discontinuousgalerkin/libdiscontinuousgalerkin.la \
		     lagrange/liblagrange.la

include $(top_srcdir)/am/global-rules
