rungekuttadir =  $(includedir)/dune/fem/solver/rungekutta
rungekutta_HEADERS = \
	basicimplicit.hh \
	basicrow.hh \
	butchertable.hh \
	explicit.hh \
	implicit.hh \
	row.hh \
	semiimplicit.hh \
	timestepcontrol.hh

noinst_LTLIBRARIES = librungekutta.la

librungekutta_la_SOURCES = butchertable.cc

include $(top_srcdir)/am/global-rules
