spacedir = $(includedir)/dune/fem/space
space_HEADERS =  \
	combinedspace.hh \
	discontinuousgalerkin.hh \
	finitevolume.hh \
	fourier.hh \
	lagrange.hh \
	padaptivespace.hh \
	rannacherturek.hh

SUBDIRS = \
	basisfunctionset \
	combinedspace \
	common \
	discontinuousgalerkin \
	finitevolume \
	fourier \
	lagrange \
	mapper \
	padaptivespace \
	rannacherturek \
	shapefunctionset \
	test

noinst_LTLIBRARIES = libspace.la

libspace_la_SOURCES =
libspace_la_LIBADD = \
	lagrange/liblagrange.la \
	shapefunctionset/libshapefunctionset.la

include $(top_srcdir)/am/global-rules
