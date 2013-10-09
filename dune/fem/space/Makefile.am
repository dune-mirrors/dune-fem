spacedir = $(includedir)/dune/fem/space

space_HEADERS = combineddiscretefunctionspace.hh \
		combinedspace.hh \
		dgspace.hh \
		discontinuousgalerkin.hh \
		finitevolume.hh \
		fourier.hh \
		fvspace.hh \
		lagrange.hh \
		lagrangespace.hh \
		padaptivespace.hh \
		rannacherturek.hh

SUBDIRS = basisfunctionset \
	  combinedspace \
	  common \
	  discontinuousgalerkin \
	  dofmapper \
	  finitevolume \
	  fourier \
	  lagrange \
	  mapper \
	  padaptivespace \
	  rannacherturek \
	  shapefunctionset \
	  test

noinst_LTLIBRARIES = libspace.la

# construct only from sub-libraries
libspace_la_SOURCES =
libspace_la_LIBADD = discontinuousgalerkin/libdiscontinuousgalerkin.la \
		     lagrange/liblagrange.la

include $(top_srcdir)/am/global-rules
